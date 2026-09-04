#!/usr/bin/env python
import argparse
import logging
import pathlib

import anndata
import geopandas
import h5py
import pandas as pd
import scipy.sparse
import spatialdata.models
import spatialdata_io

logger = logging.getLogger("baysor_denovo_to_sdata")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)


def read_table(baysor_out):
    """feature_matrix.h5 (genes x cells, CSC) plus cells.parquet -> AnnData."""
    with h5py.File(baysor_out / "feature_matrix.h5", "r") as f:
        m = f["matrix"]
        X = scipy.sparse.csc_matrix(
            (m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:])
        )
        genes = m["features/name"][:].astype(str)
        cells = m["barcodes"][:].astype(str)

    obs = pd.read_parquet(baysor_out / "cells.parquet").set_index("cell")
    obs = obs.loc[cells]
    obs["cell_id"] = obs.index.astype(str)
    obs.index = obs.index.rename(None)

    return anndata.AnnData(
        X=X.T.tocsr(), obs=obs, var=pd.DataFrame(index=pd.Index(genes, name=None))
    )


def read_boundaries(path):
    """GeoParquet polygons -> GeoDataFrame indexed by cell_id.

    Baysor emits self-intersecting rings for sparse z layers (~11% of them),
    which break unions and area computations downstream. Repaired with buffer(0);
    a few collapse to zero area and are dropped.
    """
    gdf = geopandas.read_parquet(path).rename(columns={"cell": "cell_id"})
    invalid = ~gdf.geometry.is_valid
    if invalid.any():
        gdf.loc[invalid, "geometry"] = gdf.loc[invalid, "geometry"].buffer(0)
    empty = gdf.geometry.is_empty
    logger.info(
        f"{path.name}: {len(gdf)} polygons, {invalid.sum()} repaired, "
        f"{empty.sum()} dropped as empty"
    )
    gdf = gdf[~empty]
    gdf.set_index("cell_id", drop=False, inplace=True)
    gdf.index = gdf.index.rename(None)
    return gdf


def main():
    """Convert a Baysor de novo parquet bundle to sdata.zarr."""
    parser = argparse.ArgumentParser(
        description="Convert Baysor de novo output to sdata."
    )
    parser.add_argument("data_path", help="Path to merfish output folder.")
    parser.add_argument("save_path", help="Path to Baysor_*_denovo results folder.")
    args = parser.parse_args()

    save_path = pathlib.Path(args.save_path)
    baysor_out = save_path / "baysor_out"
    assert (baysor_out / "feature_matrix.h5").exists(), "not correctly computed"

    logger.info("Loading images...")
    sdata = spatialdata_io.merscope(
        args.data_path,
        transcripts=False,
        mosaic_images=True,
        cells_boundaries=False,
    )

    logger.info("Loading boundaries and counts...")
    boundaries = read_boundaries(baysor_out / "cell_boundaries_3d.parquet")
    # layer holds the z position in microns; ZIndex is its rank among the planes.
    boundaries["ZLevel"] = boundaries["layer"].astype(float)
    levels = {v: i for i, v in enumerate(sorted(boundaries["ZLevel"].unique()))}
    boundaries["ZIndex"] = boundaries["ZLevel"].map(levels).astype(int)

    sdata["baysor_boundaries"] = spatialdata.models.ShapesModel.parse(boundaries)

    adata = read_table(baysor_out)
    adata = adata[adata.obs["cell_id"].isin(boundaries["cell_id"])].copy()
    adata.obs["region"] = pd.Categorical(["baysor_boundaries"] * adata.n_obs)
    sdata["table"] = spatialdata.models.TableModel.parse(
        adata,
        region="baysor_boundaries",
        region_key="region",
        instance_key="cell_id",
    )

    # Per z plane intensities come from intensities_3D.py, run after this.
    for i in list(sdata.images.keys()):
        del sdata[i]

    logger.info("Saving data...")
    sdata.write(str(save_path / "sdata.zarr"), overwrite=True)
    logger.info("Done.")


if __name__ == "__main__":
    main()
