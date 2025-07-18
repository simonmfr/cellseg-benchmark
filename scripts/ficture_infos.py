import logging
import os
import sys
import warnings
from os.path import exists, isdir, join
from re import split

import dask
import numpy as np
import pandas as pd
from tqdm import tqdm

dask.config.set({"dataframe.query-planning": False})

from spatialdata import SpatialData, read_zarr
from spatialdata.models import Image2DModel, ShapesModel
from spatialdata.transformations import Affine, set_transformation

from cellseg_benchmark.metrics.ficture_intensities import aggregate_channels
from cellseg_benchmark.sdata_utils import prepare_ficture

warnings.filterwarnings("ignore")

logger = logging.getLogger("ficture_infos")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

sample = sys.argv[1]
data_path = sys.argv[2]
recompute = True if sys.argv[3] == "true" else False

results_path = join(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark",
    "samples",
    sample,
    "results",
)

logger.info("Check Segementations for prior ficture results.")
compute_ficture = []
for method in os.listdir(results_path):
    if exists(join(results_path, method, "sdata.zarr")):
        if recompute:
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
stats = prepare_ficture(data_path, results_path, top_n_factors=1, logger=logger)

logger.info("Converting images")
area_covered = stats["images"] > 0
del stats
stats = prepare_ficture(data_path, results_path, logger=logger)
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
    #    del area_covered,

    logger.info("Read shapes")
    transform = pd.read_csv(
        join(data_path, "images/micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    ).values
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
                Affine(transform, input_axes=("x", "y"), output_axes=("x", "y")),
                to_coordinate_system="global",
            )
        else:
            sdata[f"boundaries_{method}"] = ShapesModel.parse(tmp[boundary_key])
    del tmp

    for key in tqdm(sdata.shapes.keys()):
        logger.info("Working on {}".format("_".join(split("_", key)[1:])))
        covered = aggregate_channels(
            sdata, image_key="ficture_image_1", shapes_key=key, mode="sum"
        )
        #        covered_weight = aggregate_channels(
        #            sdata, image_key="ficture_image_2", shapes_key=key, mode="sum"
        #        )
        #        mean = aggregate_channels(
        #            sdata, image_key="ficture_image_1", shapes_key=key, mode="average"
        #        )
        mean_weight = aggregate_channels(
            sdata, image_key="ficture_image_2", shapes_key=key, mode="average"
        )
        #        var = aggregate_channels(
        #            sdata,
        #            image_key="ficture_image_1",
        #            shapes_key=key,
        #            mode="variance",
        #            means=mean,
        #        )
        #        var_weight = aggregate_channels(
        #            sdata,
        #            image_key="ficture_image_2",
        #            shapes_key=key,
        #            mode="variance",
        #            means=mean_weight,
        #        )

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
        #        pd.DataFrame(
        #            var_weight,
        #            index=sdata[key].index.to_list(),
        #            columns=[f"fictureF21_{i}_variance_intensity_weighted" for i in range(21)],
        #        ).to_csv(
        #            join(
        #                results_path,
        #                "_".join(split("_", key)[1:]),
        #                "Ficture_stats",
        #                "vars_weight.csv",
        #            )
        #        )
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
