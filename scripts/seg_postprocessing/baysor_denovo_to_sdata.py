#!/usr/bin/env python
import argparse
import logging
import pathlib

import anndata
import geopandas
import h5py
import pandas as pd
import scipy.sparse
import sopa.aggregation
import sopa.utils
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
    which break unions and area computations downstream, so those are repaired.
    """
    gdf = geopandas.read_parquet(path).rename(columns={"cell": "cell_id"})
    invalid = ~gdf.geometry.is_valid
    if invalid.any():
        logger.info(f"repairing {invalid.sum()} invalid polygons in {path.name}")
        gdf.loc[invalid, "geometry"] = gdf.loc[invalid, "geometry"].buffer(0)
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
    boundaries_2d = read_boundaries(baysor_out / "cell_boundaries.parquet")

    sdata["baysor_boundaries"] = spatialdata.models.ShapesModel.parse(boundaries)
    sdata["baysor_boundaries_2d"] = spatialdata.models.ShapesModel.parse(boundaries_2d)

    adata = read_table(baysor_out)
    adata = adata[adata.obs["cell_id"].isin(boundaries_2d["cell_id"])].copy()
    adata.obs["region"] = pd.Categorical(["baysor_boundaries"] * adata.n_obs)
    sdata["table"] = spatialdata.models.TableModel.parse(
        adata,
        region="baysor_boundaries",
        region_key="region",
        instance_key="cell_id",
    )

    logger.info("Aggregating channel intensities...")
    sdata["table"].obsm["intensities"] = pd.DataFrame(
        sopa.aggregation.aggregate_channels(sdata, shapes_key="baysor_boundaries_2d"),
        columns=sopa.utils.validated_channel_names(
            sopa.utils.get_spatial_image(
                sdata, list(sdata.images.keys())[0], return_key=True
            )[1]
        ),
        index=sdata["baysor_boundaries_2d"].index.astype(str),
    ).loc[sdata["table"].obs["cell_id"]]

    for i in list(sdata.images.keys()):
        del sdata[i]
    del sdata["baysor_boundaries_2d"]

    logger.info("Saving data...")
    sdata.write(str(save_path / "sdata.zarr"), overwrite=True)
    logger.info("Done.")


if __name__ == "__main__":
    main()
