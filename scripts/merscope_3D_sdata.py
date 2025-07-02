import logging
import sys
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

data_path = sys.argv[1]
save_path = sys.argv[2]

assert any([".vzg" in file for file in listdir(save_path)]), "not correctly computed"
logger.info("Loading data")
sdata = merscope(
    data_path,
    transcripts=False,
    mosaic_images=False,
    cells_boundaries=False,
    vpt_outputs={
        "cell_by_gene": Path(join(save_path, "analysis_outputs", "cell_by_gene.csv")),
        "cell_metadata": Path(join(save_path, "analysis_outputs", "cell_metadata.csv")),
        "cell_boundaries": Path(
            join(save_path, "analysis_outputs", "cellpose2_micron_space.parquet")
        ),
    },
    z_layers=[0, 1, 2, 3, 4, 5, 6],
)
boundaries = read_parquet(
    join(save_path, "analysis_outputs", "cellpose2_micron_space.parquet"),
    columns=("ID", "EntityID", "ZIndex", "ZLevel", "Geometry"),
)
boundaries.rename_geometry("geometry", inplace=True)
sdata["boundaries_vpt_3D"] = ShapesModel.parse(boundaries)
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
sdata.write(join(save_path, "sdata.zarr"), overwrite=True)
logger.info("Done")
