import dask
dask.config.set({'dataframe.query-planning': False})

from cellseg_benchmark.metrics.ficture_intensities import aggregate_channels
from cellseg_benchmark.sdata_utils import prepare_ficture
from os.path import join
import numpy as np
from spatialdata import read_zarr, SpatialData
from spatialdata.models import Image2DModel
import sys
import pandas as pd
from re import split
import os
import warnings
import logging

warnings.filterwarnings("ignore")

logger = logging.getLogger("ficture_infos")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

sample = sys.argv[1]
data_path = sys.argv[2]
data_dir = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
sdata_path = join(data_dir, "samples", sample)

logger.info(f"Preparing Ficture")
stats = prepare_ficture(data_path, sdata_path, logger=logger)

logger.info(f"Converting images")
area_covered = (stats[0]>0)
stats[0] = (stats[0].astype(np.float16)/(np.finfo(np.float16).max.astype(np.uint16) - 5)) #reverse ficture normalization
logger.info("saving general stats")
pd.DataFrame(
    data=np.concat([area_covered.sum(axis=(1,2))[np.newaxis, :], stats[0].sum(axis=(1,2))[np.newaxis, :]], axis=0),
    columns=list(range(21)), index=["factor_area", "factor_area_weighted"]
).to_csv(join(sdata_path, "results", "Ficture", "general_stats.csv"))

logger.info(f"Read shapes")
master_sdata = read_zarr(join(sdata_path, "sdata_z3.zarr"), selection=("shapes",))

logger.info(f"Build temporary SpatialData")
sdata = SpatialData()
sdata["ficture_image_1"] = Image2DModel.parse(area_covered)
sdata["ficture_image_2"] = Image2DModel.parse(stats[0])
del area_covered, stats
for key in master_sdata.shapes.keys():
    sdata[key] = master_sdata[key]
del master_sdata

for key in sdata.shapes.keys():
    logger.info("Working on {}".format(key))
    covered = aggregate_channels(sdata, image_key="ficture_image_1", shapes_key=key, mode="sum")
    covered_weight = aggregate_channels(sdata, image_key="ficture_image_2", shapes_key=key, mode="sum")
    mean = aggregate_channels(sdata, image_key="ficture_image_1", shapes_key=key, mode="average")
    mean_weight = aggregate_channels(sdata, image_key="ficture_image_2", shapes_key=key, mode="average")
    var = aggregate_channels(sdata, image_key="ficture_image_1", shapes_key=key, mode="variance", means=mean)
    var_weight = aggregate_channels(sdata, image_key="ficture_image_2", shapes_key=key, mode="variance", means=mean_weight)

    os.makedirs(join(sdata_path, "results", "_".join(split("_", key)[1:]), "Ficture_stats"), exist_ok=True)
    pd.DataFrame(
        mean,
        index=sdata[key].index.to_list(),
        columns=[f"fictureF21_{i}_mean_intensity" for i in range(21)],
    ).to_csv(join(sdata_path, "results", "_".join(split("_", key)[1:]), "Ficture_stats", "means.csv"))
    pd.DataFrame(
        mean_weight,
        index=sdata[key].index.to_list(),
        columns=[f"fictureF21_{i}_mean_intensity_weighted" for i in range(21)],
    ).to_csv(join(sdata_path, "results", "_".join(split("_", key)[1:]), "Ficture_stats", "means_weight.csv"))
    pd.DataFrame(
        var,
        index=sdata[key].index.to_list(),
        columns=[f"fictureF21_{i}_variance_intensity" for i in range(21)],
    ).to_csv(join(sdata_path, "results", "_".join(split("_", key)[1:]), "Ficture_stats", "vars.csv"))
    pd.DataFrame(
        var_weight,
        index=sdata[key].index.to_list(),
        columns=[f"fictureF21_{i}_variance_intensity_weighted" for i in range(21)],
    ).to_csv(join(sdata_path, "results", "_".join(split("_", key)[1:]), "Ficture_stats", "vars_weight.csv"))
    pd.DataFrame(
        covered,
        index=sdata[key].index.to_list(),
        columns=[f"fictureF21_{i}_area" for i in range(21)],
    ).to_csv(join(sdata_path, "results", "_".join(split("_", key)[1:]), "Ficture_stats", "area.csv"))
    pd.DataFrame(
        covered_weight,
        index=sdata[key].index.to_list(),
        columns=[f"fictureF21_{i}_area_weighted" for i in range(21)],
    ).to_csv(join(sdata_path, "results", "_".join(split("_", key)[1:]), "Ficture_stats", "area_weight.csv"))