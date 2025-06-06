import logging
import os
import sys
import warnings

from spatialdata import read_zarr

from cellseg_benchmark.adata_utils import (
    dimensionality_reduction,
    filter_cells,
    filter_genes,
    integration_harmony,
    merge_adatas,
    normalize,
)

warnings.filterwarnings("ignore")

logger = logging.getLogger("integrate_adatas")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

path = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
name = sys.argv[1]

sdata_list = []
available_names = set()
logger.info("Loading data")
for f in os.listdir(os.path.join(path, "samples")):
    if f.startswith("foxf2"):
        if f not in ["foxf2_s2_r0", "foxf2_s3_r0", "foxf2_s3_r1"]:
            sdata = read_zarr(
                os.path.join(path, "samples", f, "sdata_z3.zarr"), selection=("tables",)
            )
            current_names = list(sdata.tables.keys())
            current_names = ["_".join(name.split("_")[1:]) for name in current_names]
            available_names.update(set(current_names))
            sdata_list.append((f, sdata))

save_path = os.path.join(path, "analysis", name, "plots")
os.makedirs(save_path, exist_ok=True)

adata = merge_adatas(
    sdata_list, key=name, logger=logger, do_qc=True, save_path=save_path
)
adata = filter_cells(adata, save_path=save_path, logger=logger)
adata = filter_genes(adata, save_path=save_path, logger=logger)
adata = normalize(adata, save_path=save_path, logger=logger)
dimensionality_reduction(adata, save_path=save_path, logger=logger)
adata = integration_harmony(
    adata, batch_key="sample", save_path=save_path, logger=logger
)

os.makedirs(os.path.join(path, "analysis", name, "adatas"), exist_ok=True)
if "fov" in adata.obs.columns:
    adata.obs["fov"] = adata.obs["fov"].astype(str)
adata.write(
    os.path.join(path, "analysis", name, "adatas", "adata_integrated.h5ad.gz"),
    compression="gzip",
)
