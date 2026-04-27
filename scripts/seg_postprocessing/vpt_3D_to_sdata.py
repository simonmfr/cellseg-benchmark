import argparse
import logging
import os
import pathlib
import geopandas
import pandas as pd
import spatialdata.models
import spatialdata_io
import sopa.aggregation
import sopa.utils

logger = logging.getLogger("vpt_3D_to_sdata")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

def main():
    parser = argparse.ArgumentParser(
        description="Convert vpt 3D output to sdata."
    )
    parser.add_argument("data_path", help="Path to vpt folder.")
    parser.add_argument("save_path", help="Path to output folder.")
    args = parser.parse_args()

    save_path = pathlib.Path(args.save_path)
    data_path = pathlib.Path(args.data_path)

    assert any(
        "cell_by_gene.csv" in file
        for file in os.listdir(save_path / "analysis_outputs")
    ), "not correctly computed"

    logger.info("Loading data..")
    sdata = spatialdata_io.merscope(
        data_path,
        transcripts=False,
        mosaic_images=True,
        cells_boundaries=False,
        vpt_outputs={
            "cell_by_gene": save_path / "analysis_outputs" / "cell_by_gene.csv",
            "cell_metadata": save_path / "analysis_outputs" / "cell_metadata.csv",
            "cell_boundaries": save_path / "analysis_outputs" / "cellpose2_micron_space.parquet",
        }
    )

    sdata["table"].obs.rename(columns={"EntityID": "cell_id"}, inplace=True)

    boundaries = geopandas.read_parquet(
        save_path / "analysis_outputs" / "cellpose2_micron_space.parquet",
        columns=("ID", "EntityID", "ZIndex", "ZLevel", "Geometry"),
    )
    boundaries.rename_geometry("geometry", inplace=True)
    boundaries.rename(columns={"EntityID": "cell_id"}, inplace=True)
    boundaries.set_index("cell_id", drop=False, inplace=True)
    boundaries.index = boundaries.index.rename(None)

    boundaries_2 = boundaries.copy()
    boundaries_2 = boundaries_2[["cell_id", "geometry"]]
    boundaries_2 = boundaries_2.dissolve(by="cell_id")
    boundaries_2.index = boundaries_2.index.rename(None)

    sdata["boundaries_vpt_3D"] = spatialdata.models.ShapesModel.parse(boundaries)
    sdata["boundaries_vpt_2D"] = spatialdata.models.ShapesModel.parse(boundaries_2)

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

    sdata["table"].obsm['intensities'] = pd.DataFrame(
        sopa.aggregation.aggregate_channels(sdata, shapes_key="boundaries_vpt_2D"),
        columns=sopa.utils.validated_channel_names(
            sopa.utils.get_spatial_image(sdata, list(sdata.images.keys())[0], return_key=True)[1]
        ),
        index=sdata["boundaries_vpt_2D"].index.astype(str)
    )

    for i in list(sdata.images.keys()):
        del sdata[i]
    del sdata["boundaries_vpt_2D"]

    logger.info("Saving data...")
    sdata.write(str(save_path / "sdata.zarr"), overwrite=True)
    logger.info("Done.")

if __name__ == "__main__":
    main()