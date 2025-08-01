import argparse
import logging
import os
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

from cellseg_benchmark._constants import (
    cell_type_colors,
    contamination_markers,
    selected_EC_subtypes,
)
from cellseg_benchmark.adata_utils import (
    clean_pca_umap,
    dimensionality_reduction,
    integration_harmony,
    normalize_counts,
)
from cellseg_benchmark.cell_annotation_utils import (
    annotate_cells_by_score,
    flag_contamination,
    score_cell_types,
)

today = date.today().strftime("%Y%m%d")

# warnings.filterwarnings("ignore")

# Logger setup
logger = logging.getLogger("vascular_subclustering")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

# CLI args
parser = argparse.ArgumentParser(
    description="Identify and annnotate vascular subtypes, including EC zonation. Requires prio (pan-)EC annotation in .obs of adata_integrated.h5ad"
)
parser.add_argument("cohort", help="Cohort name, e.g., 'foxf2'")
parser.add_argument(
    "seg_method", help="Segmentation method, e.g., 'Cellpose_1_nuclei_model'"
)
args = parser.parse_args()

# Paths
base_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
method_path = base_path / "analysis" / args.cohort / args.seg_method
method_path.mkdir(parents=True, exist_ok=True)
plot_path = method_path / "plots" / "vascular-subtypes"
plot_path.mkdir(parents=True, exist_ok=True)
plot_path_ecs = method_path / "plots" / "ec-zonation"
plot_path_ecs.mkdir(parents=True, exist_ok=True)
results_path = method_path / "vascular-subtypes"
results_path.mkdir(parents=True, exist_ok=True)
plot_path_spat = plot_path / "Spatial"
plot_path_spat.mkdir(parents=True, exist_ok=True)

# Load merged adata file (integrated samples)
logger.info("Loading integrated adata...")
adata = sc.read_h5ad(os.path.join(method_path, "adatas", "adata_integrated.h5ad.gz"))

cell_type_column = "cell_type_mmc_raw_revised"

# Use the dictionary keys as the cell type order
adata.obs[cell_type_column] = pd.Categorical(
    adata.obs[cell_type_column], categories=list(cell_type_colors.keys())
)
adata.obs[cell_type_column] = adata.obs[cell_type_column].cat.remove_unused_categories()
# set colors for allen_cell_type
colors = [cell_type_colors[cat] for cat in adata.obs[cell_type_column].cat.categories]
adata.uns[cell_type_column + "_colors"] = colors

logger.info("Subsetting vascular cells...")
sub_adata = adata[
    adata.obs[cell_type_column].isin(["ECs", "SMCs", "Pericytes", "VLMCs", "BAMs"])
].copy()
del adata
with rc_context({"figure.figsize": (8, 7)}):
    sc.pl.embedding(
        sub_adata, basis="X_umap_harmony_20_50", color=cell_type_column, show=False
    )
    plt.savefig(
        os.path.join(plot_path, "UMAP_integrated_harmony_subset.png"),
        dpi=200,
        bbox_inches="tight",
    )
    plt.close()

logger.info("Filtering out contaminated cells...")
# Additional filtering for cotaminated cells is applied here. This is to ensure high quality vascular subclusters.
sub_adata = flag_contamination(
    sub_adata,
    contamination_markers,
    layer="volume_log1p_norm",
    absolute_min=1,
    z_threshold=2,
    logger=logger,
)

with rc_context({"figure.figsize": (6, 6)}):
    sc.pl.embedding(
        sub_adata,
        wspace=0.2,
        ncols=5,
        basis="X_umap_harmony_20_50",
        color=[
            "contaminated",
            "contaminated_Neuronal",
            "contaminated_Oligodendrocytes",
            "contaminated_Astrocytes",
            "contaminated_Immune",
        ],
        size=50000 / sub_adata.shape[0],
        show=False,
    )
    plt.savefig(
        os.path.join(plot_path, "UMAP_cell_type_contaminations.png"),
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()

sub_adata = sub_adata[~sub_adata.obs["contaminated"]].copy()

logger.info("Re-processing vascular subset...")
batch_key = "sample"
adata = sub_adata.copy()
adata = clean_pca_umap(adata, logger=logger)
del sub_adata

adata.X = adata.layers["counts"]
assert np.issubdtype(adata.X.dtype, np.integer)

adata = normalize_counts(
    adata, save_path=plot_path, seg_method=args.seg_method, logger=logger
)
adata = dimensionality_reduction(
    adata, save_path=plot_path, point_size_factor=130000, logger=logger
)
adata = integration_harmony(
    adata,
    batch_key="sample",
    n_neighbors=10,
    n_pcs=20,
    save_path=plot_path,
    point_size_3d=None,
    point_alpha_3d=0.05,
    logger=logger,
)

logger.info("Annotating EC zonation subtypes...")

# Step 1: score marker genes
adata = score_cell_types(
    adata,
    marker_genes_dict=selected_EC_subtypes,
    top_n_genes=10,
    layer="volume_log1p_norm",
    logger=None,
)

# Step 2: assign each cell by max scoring subtype
adata = annotate_cells_by_score(
    adata,
    marker_dict=selected_EC_subtypes,
    out_col="ec_zonation",
    score_threshold=0.15,
    logger=None,
)

with rc_context({"figure.figsize": (7, 6)}):
    sc.pl.embedding(
        adata,
        basis="X_umap_harmony_10_20",
        color="ec_zonation",
        size=320000 / adata.shape[0],
        show=False,
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(plot_path_ecs, "UMAP_integrated_harmony_EC_scoring.png"),
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()

counts = adata.obs["ec_zonation"].value_counts()
logger.info(f"Value counts: {counts}")

logger.info("Subsetting and reprocessing EC subtypes...")
ec_adata = adata[adata.obs["ec_zonation"].isin(["capECs", "aECs", "vECs"])].copy()
logger.info(f"Cells removed (non-ECs): {adata.shape[0] - ec_adata.shape[0]}")
logger.info(f"Remaining ECs: {ec_adata.shape[0]}")

ec_adata = clean_pca_umap(ec_adata)
ec_adata.X = ec_adata.layers["counts"]
assert np.issubdtype(ec_adata.X.dtype, np.integer)

ec_adata = normalize_counts(
    ec_adata, save_path=plot_path_ecs, seg_method=args.seg_method, logger=logger
)
ec_adata = dimensionality_reduction(
    ec_adata, save_path=plot_path_ecs, point_size_factor=130000, logger=logger
)

# Use the dictionary keys as the cell type order
ec_adata.obs["ec_zonation"] = pd.Categorical(
    ec_adata.obs["ec_zonation"], categories=list(cell_type_colors.keys())
)
ec_adata.obs["ec_zonation"] = ec_adata.obs["ec_zonation"].cat.remove_unused_categories()
# set colors for allen_cell_type
colors = [cell_type_colors[cat] for cat in ec_adata.obs["ec_zonation"].cat.categories]
ec_adata.uns["ec_zonation_colors"] = colors

with rc_context({"figure.figsize": (8, 7)}):
    sc.pl.embedding(
        ec_adata,
        basis="X_umap_10_20",
        color="ec_zonation",
        size=320000 / ec_adata.shape[0],
        show=False,
    )
    plt.savefig(
        os.path.join(plot_path_ecs, "UMAP_reprocessed_ECs.png"),
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()

# 3d umap
sc.pp.neighbors(ec_adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(
    ec_adata, neighbors_key="neighbors", key_added="X_umap_10_20_3D", n_components=3
)
with rc_context({"figure.figsize": (8, 8)}):
    sc.pl.embedding(
        ec_adata,
        basis="X_umap_10_20_3D",
        color=["sample", "ec_zonation", "cell_type_mmc_raw_revised"],
        alpha=0.2,
        wspace=0.1,
        projection="3d",
        show=False,
    )
    plt.savefig(
        os.path.join(plot_path_ecs, "UMAP_reprocessed_ECs_3D.png"),
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()

with rc_context({"figure.figsize": (7, 6)}):
    sc.pl.embedding(
        ec_adata,
        basis="X_umap_10_20",
        color=["n_counts", "n_genes", "volume"],
        size=320000 / ec_adata.shape[0],
        cmap="turbo",
        show=False,
    )
    plt.savefig(
        os.path.join(plot_path_ecs, "UMAP_reprocessed_ECs_metadata.png"),
        dpi=100,
        bbox_inches="tight",
    )
    plt.close()

df = pd.crosstab(ec_adata.obs["condition"], ec_adata.obs["ec_zonation"], margins=True)
logger.info(f"Number of EC subtypes: \n{df}")

logger.info("Merging EC zonation subtypes back into vascular dataset...")
cell_type_column_zonation = "cell_type_incl_zonation"

# transfer labels
adata.obs = adata.obs.merge(
    ec_adata.obs[["ec_zonation"]], left_index=True, right_index=True, how="left"
)
adata.obs[cell_type_column_zonation] = (
    adata.obs["ec_zonation"]
    .astype("object")
    .fillna(adata.obs[cell_type_column].astype("object"))
)

missed_ECs_mask = adata.obs["ec_zonation"].isna() & adata.obs[cell_type_column].isin(
    ["ECs"]
)
adata.obs.loc[missed_ECs_mask, cell_type_column_zonation] = "otherECs"

# convert back to categorical
adata.obs[cell_type_column_zonation] = adata.obs[cell_type_column_zonation].astype(
    "category"
)
adata.obs.drop(columns="ec_zonation", inplace=True)

logger.info("Exporting UMAPs...")

# Use the dictionary keys as the cell type order
adata.obs[cell_type_column_zonation] = pd.Categorical(
    adata.obs[cell_type_column_zonation], categories=list(cell_type_colors.keys())
)
adata.obs[cell_type_column_zonation] = adata.obs[
    cell_type_column_zonation
].cat.remove_unused_categories()

# set colors for allen_cell_type
colors = [
    cell_type_colors[cat] for cat in adata.obs[cell_type_column_zonation].cat.categories
]
adata.uns[cell_type_column_zonation + "_colors"] = colors

with rc_context({"figure.figsize": (7, 6)}):
    sc.pl.embedding(
        adata,
        basis="X_umap_harmony_10_20",
        color=cell_type_column_zonation,
        size=320000 / adata.shape[0],
        show=False,
    )
    plt.savefig(
        os.path.join(plot_path, "UMAP_integrated_harmony_incl_zonation.png"),
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()

df = pd.crosstab(
    ec_adata.obs["condition"], ec_adata.obs[cell_type_column_zonation], margins=True
)
logger.info(f"Total cell numbers: \n{df}")
logger.info("Exporting cell type proportions to csv...")

df = (
    pd.crosstab(
        [adata.obs["condition"], adata.obs["sample"]],
        adata.obs[cell_type_column_zonation],
        normalize="index",
        margins=True,
    )
    .mul(100)
    .round(1)
)
df.to_csv(os.path.join(results_path, f"{today}_fraction_vasc_cells_per_condition.csv"))

df = (
    pd.crosstab(
        adata.obs["condition"],
        adata.obs[cell_type_column_zonation],
        normalize="index",
        margins=True,
    )
    .mul(100)
    .round(1)
)
df.to_csv(
    os.path.join(results_path, f"{today}_fraction_vasc_cells_per_condition_sum.csv")
)

logger.info("Exporting spatial plots...")
for sa in adata.obs["sample"].unique():
    adata_temp = adata[adata.obs["sample"] == sa]
    with rc_context({"figure.figsize": (10, 10)}):
        plot = sc.pl.embedding(
            adata_temp,
            basis="spatial_microns",
            color=cell_type_column_zonation,
            size=9,
            title=sa,
            show=False,
        )
    plt.savefig(
        os.path.join(plot_path_spat, f"Spatial_{sa}.png"), dpi=300, bbox_inches="tight"
    )
    plt.close()
    del adata_temp

logger.info("Saving re-processed vascular subset to h5ad...")
if "fov" not in adata.obs.columns:
    adata.obs["fov"] = ""
adata.obs["fov"] = adata.obs["fov"].astype(str)
output_path = method_path / "adatas"
output_path.mkdir(parents=True, exist_ok=True)
adata.write(output_path / "adata_integrated_vascular.h5ad.gz", compression="gzip")

logger.info("Done.")
