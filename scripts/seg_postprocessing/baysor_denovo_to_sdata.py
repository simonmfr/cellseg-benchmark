#!/usr/bin/env python
import argparse
import collections
import itertools
import logging
import pathlib
import subprocess

import anndata
import geopandas
import h5py
import numpy as np
import pandas as pd
import scipy.sparse
import scipy.spatial
import shapely
import sopa.aggregation
import sopa.io.explorer
import sopa.utils
import spatialdata
import spatialdata.models
import spatialdata.transformations
import spatialdata_io

logger = logging.getLogger("baysor_denovo_to_sdata")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)


def read_boundaries(path):
    """GeoParquet polygons -> GeoDataFrame indexed by cell_id."""
    gdf = geopandas.read_parquet(path).rename(columns={"cell": "cell_id"})
    return clean_boundaries(gdf, path.name)


def clean_boundaries(gdf, name):
    """Self-intersecting rings on sparse z layers are repaired, empty ones dropped."""
    invalid = ~gdf.geometry.is_valid
    if invalid.any():
        gdf.loc[invalid, "geometry"] = gdf.loc[invalid, "geometry"].buffer(0)
    empty = gdf.geometry.is_empty
    logger.info(
        f"{name}: {len(gdf)} polygons, {invalid.sum()} repaired, "
        f"{empty.sum()} dropped as empty"
    )
    gdf = gdf[~empty]
    gdf.set_index("cell_id", drop=False, inplace=True)
    gdf.index = gdf.index.rename(None)
    return gdf


def _cell_polygon(xy, foreign):
    """Baysor cpp-0.8.3: Delaunay of the cell, minus border triangles with other molecules."""
    tris = scipy.spatial.Delaunay(xy).simplices
    edges = [
        [tuple(sorted(e)) for e in itertools.combinations(t, 2)] for t in tris.tolist()
    ]
    count = collections.Counter(e for es in edges for e in es)
    border = collections.Counter(v for e, n in count.items() if n == 1 for v in e)
    keep = np.ones(len(tris), bool)
    for _ in range(100):
        peeled = False
        for i, es in enumerate(edges):
            if (
                keep[i]
                and sum(count[e] == 1 for e in es) == 1
                and not any(
                    count[e] == 2 and border[e[0]] == border[e[1]] == 2 for e in es
                )
                and shapely.contains_xy(shapely.Polygon(xy[tris[i]]), *foreign.T).any()
            ):
                keep[i], peeled = False, True
                for e in es:
                    count[e] -= 1
                    d = {0: -1, 1: 1}.get(count[e], 0)
                    border[e[0]] += d
                    border[e[1]] += d
        if not peeled:
            break
    return shapely.coverage_union_all(shapely.polygons(xy[tris[keep]]))


def plane_boundaries(mol):
    """Baysor's 3D cell polygons on every imaged z plane instead of binned layers."""
    rows = []
    for z, m in mol.groupby("z"):
        xy = m[["x", "y"]].to_numpy()
        cell = m["cell"].where(~m["is_noise"]).to_numpy()
        order = np.argsort(xy[:, 0])
        xs = xy[order, 0]
        for c, ids in pd.Series(cell).groupby(cell).indices.items():
            if len(ids) < 3:
                continue
            lo, hi = xy[ids].min(0), xy[ids].max(0)
            box = order[
                np.searchsorted(xs, lo[0]) : np.searchsorted(xs, hi[0], "right")
            ]
            box = box[(xy[box, 1] >= lo[1]) & (xy[box, 1] <= hi[1]) & (cell[box] != c)]
            try:
                rows.append((str(c), z, _cell_polygon(xy[ids], xy[box])))
            except scipy.spatial.QhullError:
                pass
    gdf = geopandas.GeoDataFrame(rows, columns=["cell_id", "layer", "geometry"])
    return clean_boundaries(gdf, "per-plane polygons")


def main():
    """Convert a Baysor de novo parquet bundle to sdata.zarr."""
    parser = argparse.ArgumentParser(
        description="Convert Baysor de novo output to sdata."
    )
    parser.add_argument("data_path", help="Path to merfish output folder.")
    parser.add_argument("save_path", help="Path to Baysor_*_denovo results folder.")
    parser.add_argument("--explorer", action="store_true", help="Write explorer files.")
    args = parser.parse_args()

    save_path = pathlib.Path(args.save_path)
    data_path = pathlib.Path(args.data_path)
    baysor_out = save_path / "baysor_out"
    feature_matrix = baysor_out / "feature_matrix.h5"
    if not feature_matrix.exists():
        raise FileNotFoundError(
            f"Missing {feature_matrix}; Baysor output not correctly computed"
        )

    has_images = (data_path / "images").is_dir()
    if has_images:
        logger.info("Loading images...")
        sdata = spatialdata_io.merscope(
            data_path,
            transcripts=False,
            mosaic_images=True,
            cells_boundaries=False,
        )
    else:
        logger.info("No images; boundaries stay in microns, no intensities.")
        sdata = spatialdata.SpatialData()

    logger.info("Loading boundaries...")
    boundaries_3d = baysor_out / "cell_boundaries_3d.parquet"
    is_3d = boundaries_3d.exists()
    if is_3d:
        boundaries = read_boundaries(boundaries_3d)
        mol = pd.read_parquet(
            baysor_out / "segmentation.parquet",
            columns=["x", "y", "z", "cell", "is_noise"],
        )
        if mol["z"].nunique() > boundaries["layer"].nunique():
            logger.info("Baysor binned z layers; rebuilding per imaged plane...")
            boundaries = plane_boundaries(mol)
        del mol
        boundaries["ZLevel"] = boundaries["layer"].astype(float)
        boundaries["ZIndex"] = boundaries["ZLevel"].rank(method="dense").astype(int) - 1
        # Estimated from all of a cell's molecules, wider than the union of z layers.
        outlines = read_boundaries(baysor_out / "cell_boundaries.parquet")
        intensity_shapes_key = "baysor_outlines"
    else:
        boundaries = read_boundaries(baysor_out / "cell_boundaries.parquet")
        outlines = None
        intensity_shapes_key = "baysor_boundaries"

    sdata["baysor_boundaries"] = spatialdata.models.ShapesModel.parse(boundaries)
    if outlines is not None:
        sdata["baysor_outlines"] = spatialdata.models.ShapesModel.parse(outlines)

    if has_images:
        translation = pd.read_csv(
            data_path / "images" / "micron_to_mosaic_pixel_transform.csv",
            sep=" ",
            header=None,
        )
        transform = spatialdata.transformations.Affine(
            translation.to_numpy(),
            input_axes=("x", "y"),
            output_axes=("x", "y"),
        )
        spatialdata.transformations.set_transformation(
            sdata["baysor_boundaries"], transform
        )
        if outlines is not None:
            spatialdata.transformations.set_transformation(
                sdata["baysor_outlines"], transform
            )

    logger.info("Loading counts...")
    with h5py.File(baysor_out / "feature_matrix.h5", "r") as f:
        m = f["matrix"]
        counts = scipy.sparse.csc_matrix(
            (m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:])
        )
        genes = m["features/name"][:].astype(str)
        cells = m["barcodes"][:].astype(str)

    obs = pd.read_parquet(baysor_out / "cells.parquet").set_index("cell").loc[cells]
    obs["cell_id"] = obs.index.astype(str)
    obs.index = obs.index.rename(None)

    var = pd.DataFrame({"gene": genes}, index=genes)
    adata = anndata.AnnData(counts.T.tocsr(), obs=obs, var=var)
    adata = adata[adata.obs["cell_id"].isin(boundaries["cell_id"])].copy()
    adata.obs["region"] = pd.Categorical(["baysor_boundaries"] * adata.n_obs)
    centroids = sdata[intensity_shapes_key].geometry.centroid
    adata.obsm["spatial"] = (
        centroids.get_coordinates()
        .set_axis(sdata[intensity_shapes_key].index.astype(str))
        .loc[adata.obs["cell_id"]]
        .to_numpy()
    )
    sdata["table"] = spatialdata.models.TableModel.parse(
        adata,
        region="baysor_boundaries",
        region_key="region",
        instance_key="cell_id",
    )

    # per z plane intensities come from intensities_3D.py, run after this
    if has_images:
        logger.info("Aggregating channel intensities...")
        sdata["table"].obsm["intensities"] = pd.DataFrame(
            sopa.aggregation.aggregate_channels(sdata, shapes_key=intensity_shapes_key),
            columns=sopa.utils.validated_channel_names(
                sopa.utils.get_spatial_image(
                    sdata, list(sdata.images.keys())[0], return_key=True
                )[1]
            ),
            index=sdata[intensity_shapes_key].index.astype(str),
        ).reindex(sdata["table"].obs["cell_id"])

    if outlines is not None:
        del sdata["baysor_outlines"]
    if args.explorer and is_3d:
        logger.warning("Explorer output is only supported for 2D Baysor; skipping.")
    if not args.explorer or is_3d:
        for i in list(sdata.images.keys()):
            del sdata[i]

    logger.info("Saving data...")
    zarr_path = save_path / "sdata.zarr"
    sdata.write(str(zarr_path), overwrite=True)

    if args.explorer and not is_3d:
        sdata = spatialdata.read_zarr(str(zarr_path))
        sopa.io.explorer.write(
            str(save_path / "sdata.explorer"),
            sdata,
            gene_column="gene",
            ram_threshold_gb=4,
            pixel_size=1 / translation.iloc[0, 0],
        )
        subprocess.run(["rm", "-r", str(zarr_path / "images")], check=True)

    subprocess.run(["rm", "-r", str(baysor_out)], check=True)
    logger.info("Done.")


if __name__ == "__main__":
    main()
