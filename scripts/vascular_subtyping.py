import argparse
import logging
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

# Logger Setup
logger = logging.getLogger("vascular_subclustering")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

# CLI arguments
parser = argparse.ArgumentParser(
    description="Identify vascular subtypes, including EC zonation"
)
parser.add_argument("cohort", help="Cohort name, e.g., 'foxf2'")
parser.add_argument(
    "seg_method", help="Segmentation method, e.g., 'Cellpose_1_nuclei_model'"
)
args = parser.parse_args()

# Paths
base_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
method_path = base_path / "analysis" / args.cohort / args.seg_method
plot_path = method_path / "plots" / "vascular-subtypes"
plot_path_ecs = method_path / "plots" / "ec-zonation"
plot_path_spat = plot_path / "Spatial"
results_path = method_path / "vascular-subtypes"
adata_path = method_path / "adatas"

for p in [plot_path, plot_path_ecs, plot_path_spat, results_path, adata_path]:
    p.mkdir(parents=True, exist_ok=True)

# Load data
logger.info("Loading integrated adata...")
adata = sc.read_h5ad(adata_path / "adata_integrated.h5ad.gz")

# Prepare metadata
cell_type_col = "cell_type_mmc_raw_revised"
adata.obs[cell_type_col] = pd.Categorical(
    adata.obs[cell_type_col], categories=list(cell_type_colors.keys())
).remove_unused_categories()
adata.uns[f"{cell_type_col}_colors"] = [
    cell_type_colors[ct] for ct in adata.obs[cell_type_col].cat.categories
]

# Subset vascular cells
logger.info("Subsetting vascular cells...")
vascular_types = ["ECs", "SMCs", "Pericytes", "VLMCs", "BAMs"]
sub_adata = adata[adata.obs[cell_type_col].isin(vascular_types)].copy()
del adata

with rc_context({"figure.figsize": (8, 7)}):
    sc.pl.embedding(
        sub_adata, basis="X_umap_harmony_20_50", color=cell_type_col, show=False
    )
    plt.savefig(
        plot_path / "UMAP_integrated_harmony_subset.png", dpi=200, bbox_inches="tight"
    )
    plt.close()

# Flag contamination
logger.info("Filtering out contaminated cells...")
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
        plot_path / "UMAP_cell_type_contaminations.png", dpi=150, bbox_inches="tight"
    )
    plt.close()

sub_adata = sub_adata[~sub_adata.obs["contaminated"]].copy()

# Re-process subset
logger.info("Re-processing vascular subset...")
sub_adata = clean_pca_umap(sub_adata, logger=logger)
sub_adata.X = sub_adata.layers["counts"]
assert np.issubdtype(sub_adata.X.dtype, np.integer)

sub_adata = normalize_counts(
    sub_adata, save_path=plot_path, seg_method=args.seg_method, logger=logger
)
sub_adata = dimensionality_reduction(
    sub_adata, save_path=plot_path, point_size_factor=130000, logger=logger
)
sub_adata = integration_harmony(
    sub_adata,
    batch_key="sample",
    n_neighbors=10,
    n_pcs=20,
    save_path=plot_path,
    logger=logger,
)

# Annotate EC subtypes
logger.info("Annotating EC zonation subtypes...")
sub_adata = score_cell_types(
    sub_adata, selected_EC_subtypes, top_n_genes=10, layer="volume_log1p_norm"
)
sub_adata = annotate_cells_by_score(
    sub_adata, selected_EC_subtypes, out_col="ec_zonation", score_threshold=0.15
)

with rc_context({"figure.figsize": (7, 6)}):
    sc.pl.embedding(
        sub_adata,
        basis="X_umap_harmony_10_20",
        color="ec_zonation",
        size=320000 / sub_adata.shape[0],
        show=False,
    )
    plt.tight_layout()
    plt.savefig(
        plot_path_ecs / "UMAP_integrated_harmony_EC_scoring.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()

# Subset and reprocess ECs
ec_adata = sub_adata[
    sub_adata.obs["ec_zonation"].isin(["capECs", "aECs", "vECs"])
].copy()
logger.info(f"Cells removed (non-ECs): {sub_adata.shape[0] - ec_adata.shape[0]}")
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

# Plot EC zonation
ec_adata.obs["ec_zonation"] = pd.Categorical(
    ec_adata.obs["ec_zonation"], categories=list(cell_type_colors.keys())
).remove_unused_categories()
ec_adata.uns["ec_zonation_colors"] = [
    cell_type_colors[ct] for ct in ec_adata.obs["ec_zonation"].cat.categories
]

with rc_context({"figure.figsize": (8, 7)}):
    sc.pl.embedding(
        ec_adata,
        basis="X_umap_10_20",
        color="ec_zonation",
        size=320000 / ec_adata.shape[0],
        show=False,
    )
    plt.savefig(
        plot_path_ecs / "UMAP_reprocessed_ECs.png", dpi=150, bbox_inches="tight"
    )
    plt.close()

# 3D UMAP
sc.pp.neighbors(ec_adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(ec_adata, key_added="X_umap_10_20_3D", n_components=3)

with rc_context({"figure.figsize": (8, 8)}):
    sc.pl.embedding(
        ec_adata,
        basis="X_umap_10_20_3D",
        color=["sample", "ec_zonation", cell_type_col],
        projection="3d",
        alpha=0.2,
        show=False,
    )
    plt.savefig(
        plot_path_ecs / "UMAP_reprocessed_ECs_3D.png", dpi=150, bbox_inches="tight"
    )
    plt.close()

with rc_context({"figure.figsize": (7, 6)}):
    sc.pl.embedding(
        ec_adata,
        basis="X_umap_10_20",
        color=["n_counts", "n_genes", "volume_final"],
        cmap="turbo",
        size=320000 / ec_adata.shape[0],
        show=False,
    )
    plt.savefig(
        plot_path_ecs / "UMAP_reprocessed_ECs_metadata.png",
        dpi=100,
        bbox_inches="tight",
    )
    plt.close()

# Merge zonation annotations
logger.info("Merging EC zonation subtypes back into vascular dataset...")
cell_type_zonation = "cell_type_incl_zonation"

sub_adata.obs = sub_adata.obs.merge(
    ec_adata.obs[["ec_zonation"]], left_index=True, right_index=True, how="left"
)
sub_adata.obs[cell_type_zonation] = (
    sub_adata.obs["ec_zonation"].fillna(sub_adata.obs[cell_type_col]).astype("category")
)
sub_adata.obs.loc[
    sub_adata.obs["ec_zonation"].isna() & (sub_adata.obs[cell_type_col] == "ECs"),
    cell_type_zonation,
] = "otherECs"
sub_adata.obs.drop(columns="ec_zonation", inplace=True)

sub_adata.obs[cell_type_zonation] = pd.Categorical(
    sub_adata.obs[cell_type_zonation], categories=list(cell_type_colors.keys())
).remove_unused_categories()
sub_adata.uns[f"{cell_type_zonation}_colors"] = [
    cell_type_colors[ct] for ct in sub_adata.obs[cell_type_zonation].cat.categories
]

with rc_context({"figure.figsize": (7, 6)}):
    sc.pl.embedding(
        sub_adata,
        basis="X_umap_harmony_10_20",
        color=cell_type_zonation,
        size=320000 / sub_adata.shape[0],
        show=False,
    )
    plt.savefig(
        plot_path / "UMAP_integrated_harmony_incl_zonation.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()

# Export proportions
logger.info("Exporting cell type proportions to csv...")

pd.crosstab(
    [sub_adata.obs["condition"], sub_adata.obs["sample"]],
    sub_adata.obs[cell_type_zonation],
    normalize="index",
).mul(100).round(1).to_csv(results_path / "fraction_vasc_cells_per_condition.csv")

pd.crosstab(
    sub_adata.obs["condition"], sub_adata.obs[cell_type_zonation], normalize="index"
).mul(100).round(1).to_csv(results_path / "fraction_vasc_cells_per_condition_sum.csv")

logger.info("Exporting spatial plots...")
for sample in sub_adata.obs["sample"].unique():
    with rc_context({"figure.figsize": (10, 10)}):
        sc.pl.embedding(
            sub_adata[sub_adata.obs["sample"] == sample],
            basis="spatial_microns",
            color=cell_type_zonation,
            size=9,
            title=sample,
            show=False,
        )
        plt.savefig(
            plot_path_spat / f"Spatial_{sample}.png", dpi=300, bbox_inches="tight"
        )
        plt.close()

# Save
logger.info("Saving re-processed vascular subset to h5ad...")
sub_adata.obs["fov"] = sub_adata.obs.get("fov", "").astype(str)
sub_adata.write(adata_path / "adata_integrated_vascular.h5ad.gz", compression="gzip")

logger.info("Done.")
