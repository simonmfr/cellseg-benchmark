import logging
import os
import sys
import warnings
from os.path import join, isdir, exists
from pathlib import Path

import dask
import numpy as np
import pandas as pd

dask.config.set({"dataframe.query-planning": False})

from spatialdata import SpatialData, read_zarr
from spatialdata.models import Image2DModel

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

results_path = join("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark", "samples", sample, "results")

logger.info("Check Segementations for prior ficture results.")
compute_ficture = []
for dir in os.listdir(results_path):
    if exists(join(results_path, dir, "sdata.zarr")):
        if not isdir(join(results_path, dir, "Ficture_stats")):
            compute_ficture.append(dir)

recompute_gen_stats = not exists(join(results_path, "Ficture", "general_stats.csv"))

if not (compute_ficture or recompute_gen_stats):
    logger.info("All ficture information is computed.")
    quit()

logger.info("Preparing Ficture")
stats = prepare_ficture(data_path, results_path, logger=logger)

logger.info("Converting images")
area_covered = stats[0] > 0
stats[0] = stats[0].astype(np.float16) / (
    np.finfo(np.float16).max.astype(np.uint16) - 5
)  # reverse ficture normalization
if recompute_gen_stats:
    logger.info("saving general stats")
    pd.DataFrame(
        data=np.concat(
            [
                area_covered.sum(axis=(1, 2))[np.newaxis, :],
                stats[0].sum(axis=(1, 2))[np.newaxis, :],
            ],
            axis=0,
        ),
        columns=list(range(21)),
        index=["factor_area", "factor_area_weighted"],
    ).to_csv(join(results_path, "Ficture", "general_stats.csv"))

if compute_ficture:
    logger.info("Read shapes")
    master_sdata = read_zarr((join(str(Path(results_path).parent.absolute()), "sdata_z3.zarr")), selection=("shapes",))

    logger.info("Build temporary SpatialData")
    sdata = SpatialData()
    sdata["ficture_image_1"] = Image2DModel.parse(area_covered)
    sdata["ficture_image_2"] = Image2DModel.parse(stats[0])
    del area_covered, stats
    for key in master_sdata.shapes.keys():
        sdata[key] = master_sdata[key]
    del master_sdata

    #for key in sdata.shapes.keys():
    for method in compute_ficture:
        logger.info("Working on {}".format(method))
        covered = aggregate_channels(
            sdata, image_key="ficture_image_1", shapes_key=f"boundaries_{method}", mode="sum"
        )
        covered_weight = aggregate_channels(
            sdata, image_key="ficture_image_2", shapes_key=f"boundaries_{method}", mode="sum"
        )
        mean = aggregate_channels(
            sdata, image_key="ficture_image_1", shapes_key=f"boundaries_{method}", mode="average"
        )
        mean_weight = aggregate_channels(
            sdata, image_key="ficture_image_2", shapes_key=f"boundaries_{method}", mode="average"
        )
        var = aggregate_channels(
            sdata, image_key="ficture_image_1", shapes_key=f"boundaries_{method}", mode="variance", means=mean
        )
        var_weight = aggregate_channels(
            sdata,
            image_key="ficture_image_2",
            shapes_key=f"boundaries_{method}",
            mode="variance",
            means=mean_weight,
        )

        os.makedirs(
            join(results_path, method, "Ficture_stats"),
            exist_ok=True,
        )
        pd.DataFrame(
            mean,
            index=sdata[f"boundaries_{method}"].index.to_list(),
            columns=[f"fictureF21_{i}_mean_intensity" for i in range(21)],
        ).to_csv(
            join(
                results_path,
                method,
                "Ficture_stats",
                "means.csv",
            )
        )
        pd.DataFrame(
            mean_weight,
            index=sdata[f"boundaries_{method}"].index.to_list(),
            columns=[f"fictureF21_{i}_mean_intensity_weighted" for i in range(21)],
        ).to_csv(
            join(
                results_path,
                method,
                "Ficture_stats",
                "means_weight.csv",
            )
        )
        pd.DataFrame(
            var,
            index=sdata[f"boundaries_{method}"].index.to_list(),
            columns=[f"fictureF21_{i}_variance_intensity" for i in range(21)],
        ).to_csv(
            join(
                results_path,
                method,
                "Ficture_stats",
                "vars.csv",
            )
        )
        pd.DataFrame(
            var_weight,
            index=sdata[f"boundaries_{method}"].index.to_list(),
            columns=[f"fictureF21_{i}_variance_intensity_weighted" for i in range(21)],
        ).to_csv(
            join(
                results_path,
                method,
                "Ficture_stats",
                "vars_weight.csv",
            )
        )
        pd.DataFrame(
            covered,
            index=sdata[f"boundaries_{method}"].index.to_list(),
            columns=[f"fictureF21_{i}_area" for i in range(21)],
        ).to_csv(
            join(
                results_path,
                method,
                "Ficture_stats",
                "area.csv",
            )
        )
        pd.DataFrame(
            covered_weight,
            index=sdata[f"boundaries_{method}"].index.to_list(),
            columns=[f"fictureF21_{i}_area_weighted" for i in range(21)],
        ).to_csv(
            join(
                results_path,
                method,
                "Ficture_stats",
                "area_weight.csv",
            )
        )
logger.info("Finished importing ficture infos.")