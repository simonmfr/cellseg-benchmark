import dask
dask.config.set({'dataframe.query-planning': False})

from cellseg_benchmark.metrics.ficture_intensities import *
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

logger = logging.getLogger("shape_mapping")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

sample = sys.argv[1]
data_path = sys.argv[2]
data_dir = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
sdata_path = join(data_dir, "samples", sample)

stats = prepare_ficture(data_path, sdata_path, logger=logger)

area_covered = (stats[0]>0)
area_covered_weighted = (stats[0].astype('f2')/np.iinfo(np.uint16).max) #reverse ficture normalization
data_tmp = np.concat([area_covered.sum(axis=(1,2))[np.newaxis, :], area_covered_weighted.sum(axis=(1,2))[np.newaxis, :]], axis=0)
general_stats = pd.DataFrame(data=data_tmp, columns=list(range(21)), index=["factor_area", "factor_area_weighted"])
general_stats.to_csv(join(sdata_path, "results", "Ficture", "general_stats.csv"))

master_sdata = read_zarr(join(sdata_path, "sdata_z3.zarr"), selection=("shapes",))

sdata = SpatialData()
sdata["ficture_image_1"] = Image2DModel.parse(area_covered)
sdata["ficture_image_2"] = Image2DModel.parse(area_covered_weighted)
for key in master_sdata.shapes.keys():
    sdata[key] = master_sdata[key]

for key in sdata.shapes.keys():
    covered = aggregate_channels(sdata, image_key="ficture_image_1", shapes_key=key, mode="sum")
    covered_weight = aggregate_channels(sdata, image_key="ficture_image_2", shapes_key=key, mode="sum")
    mean = aggregate_channels(sdata, image_key="ficture_image_1", shapes_key=f"boundary", mode="average")
    mean_weight = aggregate_channels(sdata, image_key="ficture_image_2", shapes_key=f"boundary", mode="average")
    var = aggregate_channels(sdata, image_key="ficture_image_1", shapes_key=f"boundary", mode="variance", means=mean)
    var_weight = aggregate_channels(sdata, image_key="ficture_image_2", shapes_key=f"boundary", mode="variance", means=mean_weight)

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