#!/usr/bin/env python
"""Joint BANKSY spatial-domain clustering across all samples of one cohort.

Two neighbourhood scales are clustered independently, because one clustering
cannot resolve both macro domains and laminae: on a 25 um raster ``k_geom=15``
is a ~55 um smoothing kernel, i.e. as wide as DG-sg itself.

* ``coarse`` -> macro domains (CTX, BS, STR, fiber tracts, ...)
* ``fine``   -> laminae (cortical layers, DG-sg, CA sp)

Clusters are written back to ``.obs`` as ``banksy_{scale}_k{k}_res{res}``.
Pick one column per scale from the plots, then label it with
``brain_regions_from_banksy.py``.
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

# Feature-mixing weight. Banksy documents 0.2 for cell typing and 0.8 for spatial
# domains; lambda = 1.0 drops own expression entirely, which washes out thin but
# transcriptionally sharp structures (DG-sg, CA sp) into their surroundings.
LAMBDA = 0.8

# k_geom = c(k_mean, k_agf) spatial neighbours defining the smoothing kernel.
SCALES = {
    "coarse": {"k_geom": [15, 30], "k_neighbors": [50, 75], "resolution": [0.2, 0.4]},
    "fine": {"k_geom": [6, 12], "k_neighbors": [15, 30], "resolution": [1.0, 1.5]},
}

R_SETUP = """
suppressPackageStartupMessages({
    library(Banksy)
    library(SpatialExperiment)
    library(SummarizedExperiment)
})

# One SpatialExperiment per sample. sample_id has to be set explicitly: without
# it every object defaults to "sample01" and runBanksyPCA(group=) cannot align
# the samples it is supposed to group by.
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


def to_r(obj):
    """Convert a pandas or numpy object to R.

    Args:
        obj: Object to convert.

    Returns:
        The corresponding R object.
    """
    with localconverter(default_converter + pandas2ri.converter + numpy2ri.converter):
        return ro.conversion.py2rpy(obj)


def to_pandas(obj):
    """Convert an R object to pandas.

    Args:
        obj: R object to convert.

    Returns:
        The corresponding pandas object.
    """
    with localconverter(default_converter + pandas2ri.converter + numpy2ri.converter):
        return ro.conversion.rpy2py(obj)


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
    "%d cells x %d genes, %d samples",
    adata.n_obs,
    adata.n_vars,
    adata.obs["sample"].nunique(),
)

ro.r(R_SETUP)

logger.info("Handing expression and coordinates to R...")
coords = pd.DataFrame(adata.obsm["spatial"], columns=["x", "y"], index=adata.obs_names)
dense = adata.layers["volume_log1p_norm"].T.toarray()
ro.globalenv["data"] = to_r(dense)
del dense
ro.globalenv["coords"] = to_r(coords)
ro.globalenv["genes"] = to_r(adata.var_names)
ro.globalenv["cells"] = to_r(adata.obs_names)
ro.globalenv["key"] = to_r(adata.obs["sample"].astype(str))
ro.r("""
rownames(data) <- genes
colnames(data) <- cells
se_list <- build_se_list(data, coords, key)
rm(data, coords, genes, cells); invisible(gc())
""")

cluster_keys = []
for scale, params in SCALES.items():
    logger.info(
        "BANKSY %s scale: k_geom=%s, lambda=%s", scale, params["k_geom"], LAMBDA
    )
    ro.globalenv["k_geom"] = ro.IntVector(params["k_geom"])
    ro.globalenv["k_neighbors"] = ro.IntVector(params["k_neighbors"])
    ro.globalenv["resolution"] = ro.FloatVector(params["resolution"])
    ro.globalenv["lambda_"] = ro.FloatVector([LAMBDA])
    ro.r("clusters <- banksy_run(se_list, k_geom, lambda_, k_neighbors, resolution)")

    # colData is ordered by sample after cbind, so join on the index, never reassign it.
    df = to_pandas(ro.globalenv["clusters"])
    df.columns = [
        re.sub(r"^clust_M\d+_lam[\d.]+_", f"banksy_{scale}_", c) for c in df.columns
    ]
    adata.obs = adata.obs.join(df)
    for col in df.columns:
        if adata.obs[col].isna().any():
            raise RuntimeError(f"{col}: cells missing a cluster after join.")
        adata.obs[col] = adata.obs[col].astype(str).astype("category")
    cluster_keys += list(df.columns)
    ro.r("rm(clusters); invisible(gc())")

logger.info("Plotting %d cluster columns...", len(cluster_keys))
for key in cluster_keys:
    plot_spatial_multiplot(
        adata, key, save_path=str(save_folder / "plots"), save_name=f"{key}.png"
    )

out = save_folder / f"spatial_reg_adata_{args.cohort}.h5ad"
adata.write(out)
logger.info("Wrote %s", out)
