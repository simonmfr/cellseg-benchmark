#!/usr/bin/env python
"""Joint BANKSY spatial-domain clustering across all samples of one cohort.

Two neighbourhood scales are clustered independently, because one kernel cannot
resolve both macro domains (CTX, BS, STR) and laminae (cortical layers, DG-sg,
CA sp). Clusters land in .obs as banksy_{scale}_k{k}_res{res}; pick one column
per scale from the plots, then label it with brain_regions_from_banksy.py.
"""

import argparse
import logging
import os
import pathlib
import re

import pandas as pd
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scanpy as sc
from rpy2.robjects import default_converter, numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter

from cellseg_benchmark.adata_utils import plot_spatial_multiplot

LAMBDA = 0.8  # 0.2 = cell typing, 0.8 = spatial domains (Banksy)

# k_geom = c(k_mean, k_agf) spatial neighbours, i.e. the smoothing kernel.
SCALES = {
    "coarse": {
        "k_geom": "c(15, 30)",
        "k_neighbors": "c(50, 75)",
        "resolution": "c(0.2, 0.4)",
    },
    "fine": {
        "k_geom": "c(6, 12)",
        "k_neighbors": "c(15, 30)",
        "resolution": "c(1.0, 1.5)",
    },
}

CONVERTER = default_converter + pandas2ri.converter + numpy2ri.converter

R_SETUP = """
suppressPackageStartupMessages({
    library(Banksy)
    library(SpatialExperiment)
    library(SummarizedExperiment)
})

# sample_id has to be set: it otherwise defaults to "sample01" for every sample,
# and runBanksyPCA(group = "sample_id") can no longer tell them apart.
build_se_list <- function(data, coords, key) {
    data <- as(data, "sparseMatrix")
    coords <- as.matrix(coords)
    stopifnot(ncol(data) == length(key))
    idx <- split(seq_along(key), key)
    se_list <- lapply(names(idx), function(nm) {
        SpatialExperiment(
            assay = list(normalized = data[, idx[[nm]], drop = FALSE]),
            spatialCoords = coords[idx[[nm]], , drop = FALSE],
            sample_id = nm
        )
    })
    names(se_list) <- names(idx)
    se_list
}

# Neighbourhood features per sample, then joint PCA and Leiden over all samples.
banksy_run <- function(se_list, k_geom, lambda, k_neighbors, resolution) {
    se <- lapply(se_list, computeBanksy, assay_name = "normalized",
                 compute_agf = TRUE, k_geom = k_geom)
    se <- do.call(cbind, se)
    invisible(gc())
    se <- runBanksyPCA(se, use_agf = TRUE, lambda = lambda,
                       group = "sample_id", seed = 1000)
    se <- clusterBanksy(se, use_agf = TRUE, lambda = lambda,
                        resolution = resolution, k_neighbors = k_neighbors,
                        seed = 1000)
    se <- connectClusters(se)
    cd <- colData(se)
    as.data.frame(cd[, grep("^clust", colnames(cd)), drop = FALSE])
}
"""

logger = logging.getLogger("banksy")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)
rcb.logger.handlers = logger.handlers
rcb.logger.setLevel(logging.ERROR)

parser = argparse.ArgumentParser(description="Joint BANKSY clustering per cohort.")
parser.add_argument("cohort", help="Cohort name, e.g. 'foxf2'.")
parser.add_argument("adata_path", help="Integrated adata of a raster segmentation.")
parser.add_argument("save_folder", help="Folder for the clustered adata and plots.")
args = parser.parse_args()

save_folder = pathlib.Path(args.save_folder)
(save_folder / "plots").mkdir(parents=True, exist_ok=True)

if "SLURM_CPUS_PER_TASK" in os.environ:
    sc.settings.n_jobs = int(os.environ["SLURM_CPUS_PER_TASK"])

logger.info("Loading %s", args.adata_path)
adata = sc.read_h5ad(args.adata_path)
logger.info(
    "%d cells, %d genes, %d samples",
    adata.n_obs,
    adata.n_vars,
    adata.obs["sample"].nunique(),
)

ro.r(R_SETUP)
logger.info("Handing expression and coordinates to R...")
with localconverter(CONVERTER):
    ro.globalenv["data"] = adata.layers["volume_log1p_norm"].T.toarray()
    ro.globalenv["coords"] = pd.DataFrame(
        adata.obsm["spatial"], columns=["x", "y"], index=adata.obs_names
    )
    ro.globalenv["genes"] = adata.var_names
    ro.globalenv["cells"] = adata.obs_names
    ro.globalenv["key"] = adata.obs["sample"].astype(str)
ro.r("""
rownames(data) <- genes
colnames(data) <- cells
se_list <- build_se_list(data, coords, key)
rm(data, coords, genes, cells); invisible(gc())
""")

cluster_keys = []
for scale, p in SCALES.items():
    logger.info("BANKSY %s: k_geom=%s, lambda=%s", scale, p["k_geom"], LAMBDA)
    with localconverter(CONVERTER):
        df = ro.r(
            f"banksy_run(se_list, {p['k_geom']}, {LAMBDA},"
            f" {p['k_neighbors']}, {p['resolution']})"
        )

    df.columns = [
        re.sub(r"^clust_M\d+_lam[\d.]+_", f"banksy_{scale}_", c) for c in df.columns
    ]
    # colData is ordered by sample after cbind, so join on the index.
    adata.obs = adata.obs.join(df)
    for col in df.columns:
        if adata.obs[col].isna().any():
            raise RuntimeError(f"{col}: cells missing a cluster after join.")
        adata.obs[col] = adata.obs[col].astype(str).astype("category")
    cluster_keys += list(df.columns)
    ro.r("invisible(gc())")

logger.info("Plotting %d cluster columns...", len(cluster_keys))
for key in cluster_keys:
    plot_spatial_multiplot(
        adata, key, save_path=str(save_folder / "plots"), save_name=f"{key}.png"
    )

out = save_folder / f"spatial_reg_adata_{args.cohort}.h5ad"
adata.write(out)
logger.info("Wrote %s", out)
