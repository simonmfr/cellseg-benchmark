import argparse
import gzip
import io
import logging
import pathlib

import geopandas as gpd
import pandas as pd
import sopa
import spatialdata as sd
import spatialdata_io

logger = logging.getLogger("intensities_3D")
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

def determine_image(sdata, idx):
    for image_name in sdata.images.keys():
        if image_name.split("_")[-1]==idx:
            return image_name
    raise ValueError("Image Level not found")

def main():
    parser = argparse.ArgumentParser(
        description=(
            "rerun image intensity computation for 3D segmentation methods and average the results per cell."
        )
    )
    parser.add_argument("sdata_path", type=str, help="Path to sdata.")
    parser.add_argument("data_path", type=str, help="Path to sis_out.")
    parser.add_argument("--method", type=str, default=None, help="method name.")
    parser.add_argument("--boundary_path", type=str, default=None, help="Path to boundaries file.")
    parser.add_argument("--boundary_key", type=str, default=None, help="Key of 3D boundary in sdata.")
    args = parser.parse_args()

    sdata_path = pathlib.Path(args.sdata_path)

    logger.info("Loading images…")
    tmp_sdata = spatialdata_io.merscope(
        args.data_path,
        transcripts=True,
        mosaic_images=True,
        cells_boundaries=False,
        z_layers=[0, 1, 2, 3, 4, 5, 6]
    )
    logger.debug(tmp_sdata)

    transform = sd.transformations.get_transformation(tmp_sdata[list(tmp_sdata.points.keys())[0]])
    logger.debug(transform)

    logger.info("Loading sdata…")
    sdata = sd.read_zarr(sdata_path / "sdata.zarr")
    logger.debug(sdata)

    logger.info("Loading boundaries…")
    logger.debug(f"boundary_key: {args.boundary_key}, method: {args.method}, boundary_path: {args.boundary_path}")
    if args.boundary_key is not None:
        boundaries = sdata[args.boundary_key]
    elif args.method is not None:
        logger.info(f"Loading boundaries from source file")
        if args.method.startswith("vpt_3D"):
            if "boundaries_vpt_3D" in sdata.shapes.keys():
                logger.debug("Loading boundaries_vpt_3D for vpt_3D by default key")
                boundaries = sdata["boundaries_vpt_3D"]
            else:
                logger.debug("Loading boundaries for vpt_3D by path")
                assert args.boundary_path.endswith(".parquet"), "vpt_3D shape files end with .parquet"
                boundaries = gpd.read_parquet(
                    args.boundary_path,
                    columns=("ID", "EntityID", "ZIndex", "ZLevel", "Geometry"),
                )
                boundaries.rename(columns={"EntityID": "cell_id"}, inplace=True)
                boundaries.set_index("cell_id", drop=False, inplace=True)
                boundaries.index = boundaries.index.rename(None)
        elif args.method.startswith("Proseg_3D"):
            logger.debug("Loading boundaries for Proseg_3D by path")
            assert args.boundary_path.endswith(".geojson.gz") or args.boundary_path.endswith(".geojson"), "This is not the Proseg 3D shapes file."
            with gzip.open(args.boundary_path, "rt", encoding="utf-8") as f:
                geojson_text = f.read()
            boundaries = gpd.read_file(io.StringIO(geojson_text))
            boundaries = boundaries.merge(sdata["table"].obs[["cell", "cell_id"]], on="cell")
        elif args.method == "Watershed_Merlin":
            if "boundaries_vpt_3D" in sdata.shapes.keys():
                logger.debug("Loading boundaries_vpt_3D for Watershed_Merlin by default key")
                boundaries = sdata["boundaries_vpt_3D"]
            else:
                logger.debug("Loading boundaries for Watershed_Merlin by path")
                boundaries = spatialdata_io.merscope(
                                args.data_path,
                                transcripts=False,
                                mosaic_images=False,
                                cells_boundaries=True
                            )
        else:
            raise NotImplementedError("Please either provide keys to the 3D boundaries in the sdata or an implemented method name.")
    else:
        raise ValueError("Please either provide keys to the sdata or a method name together with a path to the boundaries file.")

    assert sdata['table'].uns["spatialdata_attrs"]["instance_key"] in boundaries.columns, "instance_key not found in boundaries dataframe."
    assert 'ZIndex' in boundaries.columns, "ZIndex not found in boundaries dataframe."

    logger.info("Build sdata with boundaries for every level…")
    sdata_working = sd.SpatialData()
    for key, value in tmp_sdata.images.items():
        sdata_working[key] = value
    for idx in boundaries['ZIndex'].unique():
        sdata_working[f'boundaries_z{idx}'] = sd.models.ShapesModel.parse(
            boundaries[boundaries['ZIndex'] == idx],
        )
        sd.transformations.set_transformation(sdata_working[f'boundaries_z{idx}'], transform)
    sdata_working['table'] = sdata['table'].copy()
    logger.debug(sdata_working)

    logger.info("Compute intensities…")
    intensities = {}
    for idx in [x.split("_")[-1] for x in sdata_working.shapes.keys()]:
        sdata_working["table"].uns["spatialdata_attrs"]["region"] = f"boundaries_{idx}"
        sdata_working["table"].obs[sdata_working["table"].uns["spatialdata_attrs"]["region_key"]] = (
            f"boundaries_{idx}"
        )
        sdata_working["table"].obs[sdata_working["table"].uns["spatialdata_attrs"]["region_key"]] = (
            sdata_working["table"]
            .obs[sdata_working["table"].uns["spatialdata_attrs"]["region_key"]]
            .astype("category")
        )
        logger.debug(sdata_working['table'].uns["spatialdata_attrs"])
        logger.debug(sdata_working['table'].obs.head(5))
        logger.debug(sd.transformations.get_transformation(sdata_working[f"boundaries_{idx}"]))
        logger.debug(sd.transformations.get_transformation(sdata_working[determine_image(sdata_working, idx)]))

        intensities[idx] = pd.DataFrame(
            sopa.aggregation.aggregate_channels(sdata_working,
                                                shapes_key=f"boundaries_{idx}",
                                                image_key=determine_image(sdata_working, idx)
                                                ),
            columns=sopa.utils.validated_channel_names(
                sopa.utils.get_spatial_image(sdata_working, determine_image(sdata_working, idx), return_key=True)[1]
            ),
            index=sdata_working[f"boundaries_{idx}"].index.astype(str)
        )
    intensities_stacked = pd.concat(intensities.values(), keys=intensities.keys())

    logger.debug("sdata index:")
    logger.debug(sdata['table'].obs_names)
    logger.debug("intensities index")
    logger.debug(intensities_stacked.groupby(level=1).mean().index)
    sdata['table'].obsm['intensities'] = intensities_stacked.groupby(level=1).mean()

    logger.info("Write sdata with updated intensities.")
    sdata.write(sdata_path / "sdata_intens.zarr")

if __name__ == "__main__":
    main()