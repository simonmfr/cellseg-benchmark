import argparse
import logging
import pathlib
import anndata as ad
import geopandas as gpd
import pandas as pd
import spatialdata as sd
import spatialdata.models as sd_models

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]: %(message)s")
logger = logging.getLogger("sis_to_sdata")

def main(args):
    """Convert SIS output to a SpatialData zarr store.

    Matches cells via SpotTable_cell_id -> cell_label. Boundaries are 3D (multiple z-plane polygons per cell). Raises error if >5% of cells lack a boundary polygon.
    """
    sis_out = pathlib.Path(args.data_path, "sis_out")
    assert (sis_out / "cell_by_gene.h5ad").exists(), "cell_by_gene.h5ad not found"
    assert (sis_out / "cell_polygons.geojson").exists(), "cell_polygons.geojson not found"

    logger.info("Loading data...")
    adata = ad.read_h5ad(sis_out / "cell_by_gene.h5ad")
    boundaries = gpd.read_file(sis_out / "cell_polygons.geojson")

    spot_id_to_label = dict(zip(adata.obs["SpotTable_cell_id"].astype(str), adata.obs["cell_label"].astype(str)))
    boundaries["cell_label"] = boundaries["id"].astype(str).map(spot_id_to_label)
    boundaries = boundaries.dropna(subset=["cell_label", "geometry"]).copy()
    boundaries["geometry"] = boundaries.geometry.make_valid()
    boundaries.index = boundaries.pop("cell_label")

    adata.obs_names = adata.obs["cell_label"].astype(str)
    boundaries = boundaries[boundaries.index.isin(adata.obs_names)]

    missing = adata.obs_names.difference(boundaries.index)
    if len(missing) / adata.n_obs > 0.05:
        raise ValueError(f"{len(missing)}/{adata.n_obs} adata cells have no boundary (>1%)")
    if len(missing):
        logger.warning(f"{len(missing)} adata cells have no boundary, dropping from table")
        adata = adata[adata.obs_names.isin(boundaries.index)].copy()
    else:
        adata = adata.copy()

    adata.obs["region"] = pd.Categorical(["boundaries_3D"] * adata.n_obs)

    sdata = sd.SpatialData(
        shapes={"boundaries_3D": sd_models.ShapesModel.parse(boundaries)},
        tables={"table": sd_models.TableModel.parse(adata, region="boundaries_3D", region_key="region", instance_key="cell_label")},
    )

    out_path = pathlib.Path(args.data_path, "sdata.zarr")
    logger.info(f"Writing to {out_path}")
    sdata.write(out_path, overwrite=True)
    logger.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data_path", type=str, help="Path to data directory containing sis_out/")
    args = parser.parse_args()
    main(args)