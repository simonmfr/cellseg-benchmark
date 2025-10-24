# -*- coding: utf-8 -*-
# Auto-export from Jupyter notebook -> Python script with rpy2 integration
# Source notebook: banksy-setup 1 (1).ipynb
#
# Notes:
#  - Markdown cells preserved as comments.
#  - Python cells preserved.
#  - R cells (%%R / %R) are executed via rpy2: ro.r("""...""").
#  - Shell lines starting with '!' are commented to keep the script valid.
#  - If you need to pass data between Python and R, pandas2ri/numpy2ri are activated.
#
# Environment requirement:
#   conda install -c conda-forge r-base r-essentials rpy2
#
# ---------------------------------------------------------------------------
# Bootstrap rpy2
import os

import argparse
import logging
from pathlib import Path
import pandas as pd
import scanpy as sc
import squidpy as sq

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, default_converter, numpy2ri

# Activate automatic converters (pandas <-> R, numpy <-> R)
from datetime import date

today = date.today().strftime("%Y%m%d")


def R(code: str):
    """Run a multi-line R snippet safely via rpy2."""
    # Suppress warnings about UTF-8 on some systems by setting locale in R if needed.
    return ro.r(code)


def r_to_pandas(obj):
    # R object -> pandas DataFrame/Series/ndarray
    with localconverter(default_converter + pandas2ri.converter + numpy2ri.converter):
        return ro.conversion.rpy2py(obj)


def pandas_to_r(df):
    # pandas / numpy -> R object
    with localconverter(default_converter + pandas2ri.converter + numpy2ri.converter):
        return ro.conversion.py2rpy(df)


# Logger setup
logger = logging.getLogger("vascular_subclustering")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

###############################################################################
# Cell 4 - code
###############################################################################
# load R environment

rcb.logger.setLevel(logging.ERROR)
rcb.logger.handlers = logger.handlers
rcb.logger.setLevel(logger.level)

# %load_ext rpy2.ipython

###############################################################################
# Cell 5 - code
###############################################################################
R("""
suppressPackageStartupMessages({
    library(Banksy)
    library(SummarizedExperiment)
    library(SpatialExperiment)
    library(scuttle)
    library(scater)
    library(cowplot)
    library(ggplot2)
})
""")

###############################################################################
# Cell 6 - code
###############################################################################

# --- Captured outputs/comments from notebook ---
# '20250902'

###############################################################################
# Cell 7 - code
###############################################################################
# Logger setup
logger = logging.getLogger("vascular_subclustering")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

parser = argparse.ArgumentParser(
    description="DEA"
)
parser.add_argument("cohort", help="Cohort name, e.g., 'foxf2'")
parser.add_argument("--seg_method", default="Negative_Control_Rastered_25", help="Segmentation method, e.g., 'Cellpose_1_nuclei_model'")
args = parser.parse_args()

base_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
method_path = base_path / "analysis" / args.cohort / args.seg_method
data_dir = os.path.abspath("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
if "SLURM_CPUS_PER_TASK" in os.environ:
    sc.settings.n_jobs = int(os.environ["SLURM_CPUS_PER_TASK"])
    print(sc.settings.n_jobs)
logger.info("Loading integrated adata...")
adata = sc.read_h5ad(os.path.join(base_path, "analysis", args.seg_method, "adata_integrated.h5ad.gz"))

point_size_factor = 320000
celltype_col = "cell_type_mmc_raw_revised"

for cond, df in adata.obs.groupby("condition", observed=False):
    samples = df["sample"].unique()
    print(f"{cond}: {len(samples)} samples → {', '.join(samples)}")

sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)
resolutions = [0.25]
for res in resolutions:
    sc.tl.leiden(adata, key_added=f"leiden_res{res}".replace(".", "_"), resolution=res)

R("""
suppressPackageStartupMessages({
    library(Banksy)
    library(SummarizedExperiment)
    library(SpatialExperiment)
    library(scuttle)
    library(scater)
    library(cowplot)
    library(ggplot2)
})
""")

key_dict = adata.obs["sample"]
cells = adata.obs_names
genes = adata.var_names
data = adata.layers[
    "volume_log1p_norm"].T.toarray()  # see https://prabhakarlab.github.io/Banksy/articles/multi-sample.html
# get coords
coords = pd.DataFrame(adata.obsm["spatial"])
coords.columns = ["x", "y"]
coords.index = adata.obs.index
ro.globalenv['coords'] = pandas_to_r(coords)
ro.globalenv['genes'] = pandas_to_r(genes)
ro.globalenv['cells'] = pandas_to_r(cells)
ro.globalenv['data'] = pandas_to_r(data)
ro.globalenv['key_dict'] = pandas_to_r(key_dict)

R("""
# move to R
coords <- as.matrix(coords)
rownames(data) = genes
colnames(data) = cells
data <- as(data, "sparseMatrix")
# subset data matrix by sample
stopifnot(dim(data)[2] == length(key_dict))
index_list <- split(seq_along(key_dict), key_dict)
data_subset_list <- lapply(index_list, function(cols) data[, cols, drop = FALSE])
# subset coords by sample
index_list <- split(rownames(coords), key_dict)
coords_list <- lapply(index_list, function(rows) coords[rows, , drop = FALSE])
# Create a list of SpatialExperiment objects per sample
se_list <- lapply(names(data_subset_list), function(name) {
  SpatialExperiment(assay = list(normalized = data_subset_list[[name]]), 
                    spatialCoords = coords_list[[name]])
})
# Name list entries
names(se_list) <- names(data_subset_list)

rm(data)
rm(data_subset_list)

# First, compute BANKSY neighborhood feature matrices for each sample separately 
compute_agf <- TRUE
aname <- "normalized"
k_geom <- c(15, 30) # selecting k_geom https://github.com/prabhakarlab/Banksy/issues/35; used c(15, 30) in clustering individual samples
se_list <- lapply(se_list, computeBanksy, assay_name = aname, 
                  compute_agf = compute_agf, k_geom = k_geom)

# next, merge the samples to perform joint dimensional reduction and clustering
se_joint <- do.call(cbind, se_list)
rm(se_list)
invisible(gc())
lambda <- c(1.0) # can specify for multiple runs
use_agf <- c(FALSE, TRUE)
se_joint <- runBanksyPCA(se_joint, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000)
# Run UMAP on the BANKSY embedding
se_joint <- runBanksyUMAP(se_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
# obtain cluster labels for spots across all samples
res <- c(0.2, 0.4)
k_neigh <- c(30, 50, 60, 75)
se_joint <- clusterBanksy(se_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000, k_neighbors = k_neigh)
se_joint <- connectClusters(se_joint)
cnames <- colnames(colData(se_joint))
cnames <- cnames[grep("^clust", cnames)]
banksy_clusters <- colData(se_joint)[cnames]
new_cnames <- gsub("^clust_", "banksy_joined_", cnames)
new_cnames <- gsub("_lam0_", "_nonspatial_lam0_", new_cnames)
new_cnames <- gsub("_lam0.2_", "_cell_typing_lam0.2_", new_cnames)
new_cnames <- gsub("_lam0.8_", "_spatial_domains_lam0.8_", new_cnames)
colnames(banksy_clusters) <- new_cnames
""")

r_obj = ro.globalenv['banksy_clusters']
as_df = ro.r['as.data.frame']
r_df = as_df(r_obj)

# convert to pandas DataFrame
with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
    pd_df = ro.conversion.rpy2py(r_df)

adata.obs = adata.obs.join(pd_df)
for s in adata.obs['sample'].unique():
    sq.pl.spatial_scatter(
        adata[adata.obs['sample'] == s], shape=None, color=pd_df, size=0.5, library_id="spatial", figsize=(7, 7),
        wspace=0.25,
        save=f"/dss/dsshome1/00/ra87rib/cellseg-benchmark/misc/banksy_tests/{args.cohort}/banksy_align_{s}_1.png"
    )
adata.write(
    f"/dss/dsshome1/00/ra87rib/cellseg-benchmark/analysis/{args.cohort}/{args.seg_method}/spatial_reg_adata.h5ad")
