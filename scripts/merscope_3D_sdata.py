import argparse
import logging
from os import listdir
from os.path import join
from pathlib import Path

from geopandas import read_parquet
from spatialdata.models import ShapesModel
from spatialdata_io import merscope

logger = logging.getLogger("Merscope_3D")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

parser = argparse.ArgumentParser(
    description="Compute segmentation based on vpt 3D pipeline."
)
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument("save_path", help="Path to output folder.")
args = parser.parse_args()

assert any([".vzg" in file for file in listdir(args.save_path)]), (
    "not correctly computed"
)
logger.info("Loading data")
sdata = merscope(
    args.data_path,
    transcripts=False,
    mosaic_images=False,
    cells_boundaries=False,
    vpt_outputs={
        "cell_by_gene": Path(
            join(args.save_path, "analysis_outputs", "cell_by_gene.csv")
        ),
        "cell_metadata": Path(
            join(args.save_path, "analysis_outputs", "cell_metadata.csv")
        ),
        "cell_boundaries": Path(
            join(args.save_path, "analysis_outputs", "cellpose2_micron_space.parquet")
        ),
    },
    z_layers=[0, 1, 2, 3, 4, 5, 6],
)
sdata["table"].obs.rename(columns={"EntityID": "cell_id"}, inplace=True)
boundaries = read_parquet(
    join(args.save_path, "analysis_outputs", "cellpose2_micron_space.parquet"),
    columns=("ID", "EntityID", "ZIndex", "ZLevel", "Geometry"),
)
boundaries.rename_geometry("geometry", inplace=True)
boundaries.rename(columns={"EntityID": "cell_id"}, inplace=True)
boundaries.index = boundaries["cell_id"]
sdata["boundaries_vpt_3D"] = ShapesModel.parse(boundaries)
sdata["table"].uns["spatialdata_attrs"]["instance_key"] = "cell_id"
sdata["table"].uns["spatialdata_attrs"]["region"] = "boundaries_vpt_3D"
sdata["table"].obs[sdata["table"].uns["spatialdata_attrs"]["region_key"]] = (
    "boundaries_vpt_3D"
)
sdata["table"].obs[sdata["table"].uns["spatialdata_attrs"]["region_key"]] = (
    sdata["table"]
    .obs[sdata["table"].uns["spatialdata_attrs"]["region_key"]]
    .astype("category")
)

logger.info("Saving data")
sdata.write(join(args.save_path, "sdata.zarr"), overwrite=True)
logger.info("Done")
