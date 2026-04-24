import argparse
import logging
import os
import warnings
from os.path import exists, join

import dask
import numpy as np
import pandas as pd
from tqdm import tqdm

from cellseg_benchmark.metrics.ovrl import (
    compute_mean_vsi_per_polygon,
    compute_ovrl,
    plot_vsi_overview,
)

dask.config.set({"dataframe.query-planning": False})

from spatialdata import SpatialData, read_zarr, transform

from cellseg_benchmark.sdata_utils import get_2D_boundaries

warnings.filterwarnings("ignore")

logger = logging.getLogger("ficture_infos")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

parser = argparse.ArgumentParser(description="Compute Ovrlpy statistics.")
parser.add_argument("sample", help="sample name.")
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument(
    "--recompute", default=False, action="store_true", help="Recompute all stats."
)
args = parser.parse_args()

sample_path = join(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark", "samples", args.sample
)
results_path = join(
    sample_path,
    "results",
)

logger.info("Check Segementations for prior ovrlpy results.")
compute_ovrlpy = []
for method in os.listdir(results_path):
    if exists(join(results_path, method, "sdata.zarr")):
        if args.recompute:
            compute_ovrlpy.append(method)
        elif not exists(join(results_path, method, "Ovrlpy_stats", "Ovrlpy_stats.csv")):
            compute_ovrlpy.append(method)

if not compute_ovrlpy:
    logger.info("All ovrlpy information is computed.")
    quit()
else:
    logger.info(
        f"Identified {len(compute_ovrlpy)} Segmentations to compute Ovrlpy for."
    )

logger.info("Preparing Ovrlpy statistics.")
compute_ovrl(args.sample, sample_path, args.data_path, logger=logger)

if compute_ovrlpy:
    logger.info("Read shapes")
    transformation = pd.read_csv(
        join(args.data_path, "images/micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    ).values

    integrity_matrix = np.load(
        join(sample_path, "vertical_doublets_ovrlpy_output.npz")
    )["integrity"]
    signal_matrix = np.load(join(sample_path, "vertical_doublets_ovrlpy_output.npz"))[
        "signal"
    ]
    for method in tqdm(compute_ovrlpy):
        logger.info(f"Computing Ovrlpy statistics for {method}.")
        save_path = join(results_path, method, "Ovrlpy_stats")
        os.makedirs(save_path, exist_ok=True)
        sdata = SpatialData()
        tmp = read_zarr(
            join(results_path, method, "sdata.zarr"), selection=("shapes", "tables")
        )
        boundary_key = tmp[list(tmp.tables.keys())[0]].uns["spatialdata_attrs"][
            "region"
        ]
        get_2D_boundaries(method, tmp, sdata, transformation, boundary_key)
        boundary = transform(
            sdata[f"boundaries_{method}"], to_coordinate_system="micron"
        )
        res = compute_mean_vsi_per_polygon(integrity_matrix, boundary, transformation)
        res.to_csv(join(save_path, "Ovrlpy_stats.csv"))
        plot_vsi_overview(
            integrity_map=integrity_matrix,
            signal_map=signal_matrix,
            boundaries_aligned=boundary,
            vsi_mean=res["mean_integrity"].values,
            sample_name=args.sample,
            png_path=join(save_path, "ovrlpy_overview_plot.png"),
        )

logger.info("Finished computing ovrlpy.")
