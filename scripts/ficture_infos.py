import argparse
import logging
import os
import warnings
from os.path import exists, isdir, join
from re import split

import dask
import numpy as np
import pandas as pd
from tqdm import tqdm

dask.config.set({"dataframe.query-planning": False})

from spatialdata import SpatialData, read_zarr, transform
from spatialdata.models import Image2DModel, ShapesModel
from spatialdata.transformations import Affine, set_transformation, Identity

from cellseg_benchmark.metrics.ficture_intensities import _aggregate_channels_aligned
from cellseg_benchmark.sdata_utils import prepare_ficture
from cellseg_benchmark._constants import image_based

warnings.filterwarnings("ignore")

logger = logging.getLogger("ficture_infos")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

parser = argparse.ArgumentParser(description="Compute Ficture statistics.")
parser.add_argument("sample", help="sample name.")
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument(
    "--recompute", default=False, type="store_true", help="Recompute all stats."
)
args = parser.parse_args()

results_path = join(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark",
    "samples",
    args.sample,
    "results",
)

logger.info("Check Segementations for prior ficture results.")
compute_ficture = []
for method in os.listdir(results_path):
    if exists(join(results_path, method, "sdata.zarr")):
        if args.recompute:
            compute_ficture.append(method)
        elif not isdir(join(results_path, method, "Ficture_stats")):
            compute_ficture.append(method)

recompute_gen_stats = not exists(join(results_path, "Ficture", "general_stats.csv"))

if not (compute_ficture or recompute_gen_stats):
    logger.info("All ficture information is computed.")
    quit()
else:
    logger.info(
        f"Identified {len(compute_ficture)} Segmentations to compute Ficture for."
    )

logger.info("Preparing Ficture")
stats = prepare_ficture(args.data_path, results_path, top_n_factors=1, logger=logger)

logger.info("Converting images")
area_covered = stats["images"] > 0
del stats
stats = prepare_ficture(args.data_path, results_path, logger=logger)
stats["images"] = stats["images"].astype(np.float16) / (
    np.finfo(np.float16).max.astype(np.uint16) - 5
)  # reverse ficture normalization
if recompute_gen_stats:
    logger.info("saving general stats")
    pd.DataFrame(
        data=np.concat(
            [
                area_covered.sum(axis=(1, 2))[np.newaxis, :],
                #                stats[0].sum(axis=(1, 2))[np.newaxis, :],
            ],
            axis=0,
        ),
        columns=list(range(21)),
        index=["factor_area"],
        #        index=["factor_area", "factor_area_weighted"],
    ).to_csv(join(results_path, "Ficture", "general_stats.csv"))

if compute_ficture:
    logger.info("Build temporary SpatialData")
    sdata = SpatialData()
    sdata["ficture_image_1"] = Image2DModel.parse(area_covered)
    sdata["ficture_image_2"] = Image2DModel.parse(stats["images"])
    del area_covered, stats

    logger.info("Read shapes")
    transformation = pd.read_csv(
        join(args.data_path, "images/micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    ).values

    set_transformation(
        sdata[f"ficture_image_1"],
        Affine(transformation, input_axes=("x", "y"), output_axes=("x", "y")).inverse(),
        to_coordinate_system="micron",
    )
    set_transformation(
        sdata[f"ficture_image_2"],
        Affine(transformation, input_axes=("x", "y"), output_axes=("x", "y")).inverse(),
        to_coordinate_system="micron",
    )

    for method in compute_ficture:
        tmp = read_zarr(
            join(results_path, method, "sdata.zarr"), selection=("shapes", "tables")
        )
        boundary_key = tmp[list(tmp.tables.keys())[0]].uns["spatialdata_attrs"][
            "region"
        ]
        if method.startswith("vpt_3D"):
            bound = tmp[boundary_key][["EntityID", "geometry"]].dissolve(by="EntityID")
            sdata[f"boundaries_{method}"] = ShapesModel.parse(bound)
            set_transformation(
                sdata[f"boundaries_{method}"],
                Affine(transformation, input_axes=("x", "y"), output_axes=("x", "y")),
                to_coordinate_system="global",
            )
        else:
            sdata[f"boundaries_{method}"] = ShapesModel.parse(tmp[boundary_key])
        if any([method.startswith(x) for x in image_based]):
            set_transformation(
                sdata[f"boundaries_{method}"],
                Affine(transformation, input_axes=("x", "y"), output_axes=("x", "y")).inverse(),
                to_coordinate_system="micron",
            )
        else:
            set_transformation(
                sdata[f"boundaries_{method}"],
                Identity(),
                to_coordinate_system="micron",
            )
    del tmp

    image_1 = transform(sdata["ficture_image_1"], to_coordinate_system="micron")
    image_2 = transform(sdata["ficture_image_2"], to_coordinate_system="micron")
    for key in tqdm(sdata.shapes.keys()):
        logger.info("Working on {}".format("_".join(split("_", key)[1:])))
        boundary = transform(sdata[f"boundaries_{key}"], to_coordinate_system="micron")
        covered = _aggregate_channels_aligned(
            image=image_1, geo_df=boundary, mode="sum"
        )
        #        covered_weight = aggregate_channels(
        #            sdata, image_key="ficture_image_2", shapes_key=key, mode="sum"
        #        )
        #        mean = aggregate_channels(
        #            sdata, image_key="ficture_image_1", shapes_key=key, mode="average"
        #        )
        mean_weight = _aggregate_channels_aligned(
            image=image_2, geo_df=boundary, mode="average"
        )
        #        var = aggregate_channels(
        #            sdata,
        #            image_key="ficture_image_1",
        #            shapes_key=key,
        #            mode="variance",
        #            means=mean,
        #        )
        var_weight = _aggregate_channels_aligned(
            image=image_2, geo_df=boundary,
            mode="variance",
            means=mean_weight,
        )

        os.makedirs(
            join(results_path, "_".join(split("_", key)[1:]), "Ficture_stats"),
            exist_ok=True,
        )
        #        pd.DataFrame(
        #            mean,
        #            index=sdata[key].index.to_list(),
        #            columns=[f"fictureF21_{i}_mean_intensity" for i in range(21)],
        #        ).to_csv(
        #            join(
        #                results_path,
        #                "_".join(split("_", key)[1:]),
        #                "Ficture_stats",
        #                "means.csv",
        #            )
        #        )
        pd.DataFrame(
            mean_weight,
            index=sdata[key].index.to_list(),
            columns=[f"fictureF21_{i}_mean_intensity_weighted" for i in range(21)],
        ).to_csv(
            join(
                results_path,
                "_".join(split("_", key)[1:]),
                "Ficture_stats",
                "means_weight.csv",
            )
        )
        #        pd.DataFrame(
        #            var,
        #            index=sdata[key].index.to_list(),
        #            columns=[f"fictureF21_{i}_variance_intensity" for i in range(21)],
        #        ).to_csv(
        #            join(
        #                results_path,
        #                "_".join(split("_", key)[1:]),
        #                "Ficture_stats",
        #                "vars.csv",
        #            )
        #        )
        pd.DataFrame(
            var_weight,
            index=sdata[key].index.to_list(),
            columns=[f"fictureF21_{i}_variance_intensity_weighted" for i in range(21)],
        ).to_csv(
            join(
                results_path,
                "_".join(split("_", key)[1:]),
                "Ficture_stats",
                "vars_weight.csv",
            )
        )
        pd.DataFrame(
            covered,
            index=sdata[key].index.to_list(),
            columns=[f"fictureF21_{i}_area" for i in range(21)],
        ).to_csv(
            join(
                results_path,
                "_".join(split("_", key)[1:]),
                "Ficture_stats",
                "area.csv",
            )
        )
#        pd.DataFrame(
#            covered_weight,
#            index=sdata[key].index.to_list(),
#            columns=[f"fictureF21_{i}_area_weighted" for i in range(21)],
#        ).to_csv(
#            join(
#                results_path,
#                "_".join(split("_", key)[1:]),
#                "Ficture_stats",
#                "area_weight.csv",
#            )
#        )
logger.info("Finished importing ficture infos.")
