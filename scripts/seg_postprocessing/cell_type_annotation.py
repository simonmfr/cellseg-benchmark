"""Run automated cell type annotation by combining MapMyCells (MMC) and marker genes scores.

1. Load adata from sdata.zarr
2. Run MapMyCells against ABC mouse brain reference atlas; parse per-cell labels/scores
3. QC MapMyCells output: plot distributions; mark low-correlation cells as Undefined (MAD-based per group,as suggested by Allen Institute); group fine types
4. Sensitivity: annotate "mixed" cells (runner-up probability gap) and low_quality cells (MAD rule); saved as "*_mixed" and "*_low_quality"
5. Leiden clustering; assign cluster labels by MMC majority vote (from MMC-Leiden crosstab)
6. Score marker genes (ABCAtlas via sc.tl.score) and refine cluster labels using thresholds:
   - Reassign if some cell type’s cluster-mean marker-gene score ≥ 0.5 (high_threshold)
     AND it exceeds the MMC majority label’s marker-gene score by ≥ 0.25 (delta)
   - Set to "Undefined" if all cell types’ scores are < 0.5 (low_threshold)
   - Otherwise keep the MMC majority label
7. Plot UMAP and spatial plots (mixed/low-quality, annotations)
8. Export CSVs (adata_obs.csv including cell type labels).
"""

import argparse
import json
import logging
import math
import os
import warnings
from datetime import date
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from matplotlib.pyplot import rc_context
from spatialdata import read_zarr

import cellseg_benchmark.cell_annotation_utils as anno_utils
from cellseg_benchmark._constants import cell_type_colors
from cellseg_benchmark.dea_utils import add_ensembl_id

plt.rcParams["font.family"] = (
    "Arial" if "Arial" in [f.name for f in fm.fontManager.ttflist] else "sans-serif"
)
plt.rcParams["font.weight"] = "normal"
# see https://medium.com/@daimin0514/how-to-install-the-arial-font-in-linux-9e6ac76d3d9f

today = date.today().strftime("%Y%m%d")

logger = logging.getLogger("annotation")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

if "SLURM_CPUS_PER_TASK" in os.environ:
    sc.settings.n_jobs = int(os.environ["SLURM_CPUS_PER_TASK"])
    logger.info("Using SLURM_CPUS_PER_TASK={}".format(sc.settings.n_jobs))

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(
    description="Run cell type annotation using MapMyCells reference mapping to ABC mouse brain atlas."
)
parser.add_argument("sample_name", help="Sample name")
parser.add_argument("seg_method", help="Segmentation method to annotate")
parser.add_argument("data_dir", help="Base data directory")
parser.add_argument(
    "--mad_factor",
    default=3,
    type=float,
    help="MAD factor (>0) for removing outlier annotations",
)
parser.add_argument(
    "--leiden_res", default=10.0, type=float, help="Leiden clustering resolution"
)
args = parser.parse_args()

if args.mad_factor <= 0:
    parser.error("--mad_factor must be positive")

method_path = os.path.join(
    args.data_dir, "samples", args.sample_name, "results", args.seg_method
)
annotation_path = os.path.join(method_path, "cell_type_annotation")
os.makedirs(annotation_path, exist_ok=True)
mmc_dir = os.path.join(annotation_path, "mapmycells_out")
os.makedirs(mmc_dir, exist_ok=True)

logger.info("Loading data...")
adata = read_zarr(os.path.join(method_path, "sdata.zarr"))["table"]
logger.debug(f"adata columns: {adata.obs.columns}")
adata = adata[:, ~adata.var_names.str.startswith("Blank")]  # remove blank genes
adata.var["gene"] = adata.var.index

# Fix mislabeled genes in public data from Vizgen
if "VizgenMouseBrain" in args.sample_name:
    adata.var_names = adata.var_names.str.replace(r"^ADGRF3$", "Adgrf3", regex=True)
    adata.var["gene"] = adata.var["gene"].replace("ADGRF3", "Adgrf3")
    mask = adata.var_names != "missing"
    adata = adata[:, mask].copy()

logger.info("Adding ensembl IDs...")
adata.var = add_ensembl_id(
    adata.var, species="mouse", out_col="ensmus_id", logger=logger
)

pattern = f"MapMyCells_{args.sample_name}_{args.seg_method}.json"
files = [f for f in Path(mmc_dir).glob(f"*{pattern}")]
json_path = max(files, key=os.path.getmtime) if files else None

if json_path:
    logger.info(f"Using existing MapMyCells output: {json_path.name}")
else:
    logger.info("Running MapMyCells using ABC atlas as reference...")
    anno_utils.run_mapmycells(
        adata,
        sample_name=args.sample_name,
        method_name=args.seg_method,
        annotation_path=annotation_path,
        data_dir=args.data_dir,
    )
    json_path = (
        Path(mmc_dir) / f"{today}_MapMyCells_{args.sample_name}_{args.seg_method}.json"
    )

logger.info("Processing MapMyCells output...")
with open(json_path, "rb") as src:
    json_results = json.load(src)
allen_mmc_metadata = anno_utils.process_mapmycells_output(json_results)

anno_utils.plot_metric_distributions(
    allen_mmc_metadata,
    out_path=annotation_path,
    file_name="QC_raw_metric_distributions",
)

# Plot MAD thresholds for SUBC and CLAS (per cell type)
for allen_key, figsize in [("SUBC", (45, 7)), ("CLAS", (13, 7))]:
    anno_utils.plot_mad_thresholds(
        allen_mmc_metadata,
        out_path=annotation_path,
        name=f"QC_mad_threshold_{allen_key}",
        group_column=f"allen_{allen_key}",
        value_column=f"allen_cor_{allen_key}",
        mad_factor=args.mad_factor,
        figsize=figsize,
    )

logger.info("Marking low-quality cells as 'Undefined'...")
# based on correlation MAD as suggested by Allen Institute
taxonomy_levels = ["CLAS", "SUBC", "SUPT", "CLUS"]
prefixes = ["allen", "allen_runner_up_1", "allen_runner_up_2"]
for level in taxonomy_levels:
    for prefix in prefixes:
        anno_utils.mark_low_quality_mappings(
            allen_mmc_metadata,
            target_column=prefix,
            mad_factor=args.mad_factor,
            level=level,
        )

# re-group cell types from *_SUBC
for prefix in prefixes:
    for suffix in ["SUBC", "SUBC_incl_low_quality"]:
        allen_mmc_metadata[f"{prefix}_{suffix}"] = anno_utils.group_cell_types(
            allen_mmc_metadata[f"{prefix}_{suffix}"]
        )

logger.info("Annotating 'mixed' cells based on runner-up probability...")
mixed_df = anno_utils.create_mixed_cell_types(df=allen_mmc_metadata, diff_threshold=0.5)

allen_mmc_metadata = allen_mmc_metadata.merge(
    mixed_df, left_index=True, right_index=True, how="left"
)

adata.obsm["allen_cell_type_mapping"] = allen_mmc_metadata.loc[adata.obs.index]

# if matrix contains integer-like floats, convert to int64
for name, arr in [("X", adata.X), *adata.layers.items()]:
    values = arr.data if sp.issparse(arr) else np.asarray(arr)
    if np.issubdtype(values.dtype, np.floating) and np.allclose(
        values, np.round(values), rtol=0, atol=1e-8
    ):
        casted = (
            arr.astype(np.int64)
            if sp.issparse(arr)
            else arr.astype(np.int64, copy=False)
        )
        if name == "X":
            adata.X = casted
        else:
            adata.layers[name] = casted
        logger.info(f"Converting {name} (integer-like floats) to int64.")

adata = anno_utils.process_adata(adata=adata, seg_method=args.seg_method, logger=logger)

pt_size_umap = 220000 / adata.shape[0]

# plot QC metrics (mapping correlation and probability)
with rc_context({"figure.figsize": (9, 9)}):
    # Copy relevant columns from obsm to obs temporarily
    for col in adata.obsm["allen_cell_type_mapping"].columns:
        adata.obs[col] = adata.obsm["allen_cell_type_mapping"][col]
    sc.pl.umap(
        adata,
        color=[
            "allen_cor_CLAS",
            "allen_cor_SUBC",
            "allen_cor_SUPT",
            "allen_cor_CLUS",
            "allen_prob_CLAS",
            "allen_prob_SUBC",
            "allen_prob_SUPT",
            "allen_prob_CLUS",
        ],
        size=pt_size_umap,
        legend_fontoutline=2,
        legend_fontsize=10,
        ncols=4,
        vmin=0,
        vmax=1,
        cmap="cividis",
        show=False,
    )
    plt.tight_layout()
    plt.gca().set_aspect(1)
    plt.savefig(os.path.join(annotation_path, "UMAP_mapmycells_metrics.png"))
    plt.close()
    adata.obs.drop(columns=adata.obsm["allen_cell_type_mapping"].columns, inplace=True)

# plot mixed and low-quality cells on umap and spatial plot
adata.obs["cell_type_mmc_is_low_quality"] = adata.obs[
    "cell_type_mmc_incl_low_quality"
].apply(lambda x: "undefined" if x == "Undefined" else "mapped")
palette_dict = {
    "cell_type_mmc_is_mixed": {"mixed": "red", "unique": "lightgrey"},
    "cell_type_mmc_is_low_quality": {"undefined": "blue", "mapped": "lightgrey"},
}

fig, axes = plt.subplots(
    2, 2, figsize=(18, 16), gridspec_kw={"wspace": 0.15, "hspace": 0.15}
)

keys = ["cell_type_mmc_is_mixed", "cell_type_mmc_is_low_quality"]

for i, key in enumerate(keys):
    palette = palette_dict[key]

    sc.pl.umap(
        adata,
        color=key,
        size=pt_size_umap,
        legend_fontoutline=2,
        legend_fontweight="normal",
        legend_fontsize=7,
        palette=palette,
        ax=axes[0, i],
        show=False,
    )
    axes[0, i].set_aspect("equal")

    sc.pl.embedding(
        adata,
        basis="spatial",
        color=key,
        size=150000 / adata.shape[0],
        legend_loc=None,
        palette=palette,
        ax=axes[1, i],
        show=False,
    )
    axes[1, i].set_aspect("equal")

output_path = os.path.join(
    annotation_path, "UMAP_and_Spatial_mapmycells_mixed_and_undefined.png"
)
plt.savefig(output_path, dpi=200, bbox_inches="tight")
plt.close()

logger.info(
    "Number of mixed mmc: {}".format(adata.obs["cell_type_mmc_is_mixed"].value_counts())
)
logger.debug(
    "Number of mmc including mixed cells: {}".format(
        adata.obs["cell_type_mmc_incl_mixed"].value_counts()
    )
)
logger.debug(
    "Number of mmc with mixed names: {}".format(
        adata.obs["cell_type_mmc_mixed_names"].value_counts()
    )
)
logger.info(
    "Number of low-quality cells: {}".format(
        adata.obs["cell_type_mmc_is_low_quality"].value_counts()
    )
)
logger.debug(
    "Number of cells including low-quality cells: {}".format(
        adata.obs["cell_type_mmc_incl_low_quality"].value_counts()
    )
)


leiden_col = f"leiden_res{args.leiden_res}".replace(".", "_")
adata, annotation_results, mmc_leiden_crosstab = anno_utils.revise_annotations(
    adata,
    leiden_res=args.leiden_res,
    leiden_col=leiden_col,
    cell_type_colors=cell_type_colors,
    score_high_threshold=0.5,
    score_low_threshold=0.5,
    score_delta=0.25,
    top_n_genes=50,
    ABCAtlas_marker_df_path=os.path.join(
        args.data_dir,
        "misc",
        "scRNAseq_ref_ABCAtlas_Yao2023Nature",
        "marker_genes_df",
        "20250416_cell_type_markers_top50.csv",
    ),
    logger=logger,
)

adata.obs.rename(
    columns={
        "cell_type_mmc_raw_revised": "cell_type_revised",
        "cell_type_mmc_incl_mixed_revised": "cell_type_incl_mixed_revised",
        "cell_type_mmc_incl_low_quality_revised": "cell_type_incl_low_quality_revised",
    },
    inplace=True,
)

plot_keys = [
    "cell_type_mmc_raw",
    "cell_type_mmc_raw_clusters",
    "cell_type_revised",
    "cell_type_mmc_incl_mixed",
    "cell_type_mmc_incl_mixed_clusters",
    "cell_type_incl_mixed_revised",
    "cell_type_mmc_incl_low_quality",
    "cell_type_mmc_incl_low_quality_clusters",
    "cell_type_incl_low_quality_revised",
    leiden_col,
]

n_clusters = adata.obs[leiden_col].nunique()
base_colors = (
    sns.color_palette("tab20", 20)
    + sns.color_palette("tab20b", 20)
    + sns.color_palette("tab20c", 20)
    + sns.color_palette("husl", 30)
    + sns.color_palette("Set1", 9)
    + sns.color_palette("Set2", 8)
    + sns.color_palette("Set3", 12)
)
long_palette = (base_colors * ((n_clusters // len(base_colors)) + 1))[:n_clusters]

# consistent category order for cell-type keys
ct_categories = list(cell_type_colors.keys())
for k in plot_keys:
    if k != leiden_col and k in adata.obs:
        adata.obs[k] = pd.Categorical(
            adata.obs[k], categories=ct_categories
        ).remove_unused_categories()

size_umap = pt_size_umap
size_spatial = 110000 / adata.shape[0]


def _plot_row(ax_umap, ax_umap_labels, ax_spatial, key):
    """Plotting helper."""
    is_leiden = key == leiden_col
    palette = long_palette if is_leiden else cell_type_colors

    sc.pl.umap(
        adata,
        color=key,
        size=size_umap,
        palette=palette,
        legend_loc=None,
        legend_fontsize=7,
        legend_fontweight="normal",
        legend_fontoutline=2,
        ax=ax_umap,
        show=False,
    )
    ax_umap.set_aspect("equal")

    sc.pl.umap(
        adata,
        color=key,
        size=size_umap,
        palette=palette,
        legend_loc="on data",
        legend_fontsize=9,
        legend_fontweight="normal",
        legend_fontoutline=1.5,
        ax=ax_umap_labels,
        show=False,
    )
    ax_umap_labels.set_aspect("equal")

    sc.pl.embedding(
        adata,
        basis="spatial",
        color=key,
        size=size_spatial,
        palette=palette,
        legend_loc=None if is_leiden else "right margin",
        ax=ax_spatial,
        show=False,
    )
    ax_spatial.set_aspect("equal")


n_rows = len(plot_keys)
fig, axes = plt.subplots(
    n_rows,
    3,
    figsize=(22, 7.5 * n_rows),
    gridspec_kw={"wspace": 0.05, "hspace": 0.10},
)
axes = np.atleast_2d(axes)
for i, key in enumerate(plot_keys):
    _plot_row(axes[i, 0], axes[i, 1], axes[i, 2], key)

out_png = os.path.join(annotation_path, "UMAP_and_Spatial_annotation_results.png")
plt.savefig(out_png, dpi=120, bbox_inches="tight")
plt.close()

# faceted per–cell type plot
ct_col = "cell_type_revised"
spatial = adata.obsm["spatial"]
cell_types = adata.obs[ct_col].cat.categories
n_types = len(cell_types)
n_cols = math.ceil(math.sqrt(n_types))
n_rows = math.ceil(n_types / n_cols)

fig, axs = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
axs = axs.flatten() if n_types > 1 else [axs]

# Spatial limits with margin
x, y = spatial[:, 0], spatial[:, 1]
margin_x, margin_y = 0.05 * (x.max() - x.min()), 0.05 * (y.max() - y.min())
xlim, ylim = (
    (x.min() - margin_x, x.max() + margin_x),
    (y.min() - margin_y, y.max() + margin_y),
)

handles, labels = [], []
for i, ct in enumerate(cell_types):
    if i >= len(axs):
        break
    mask = adata.obs[ct_col] == ct
    if not mask.any():
        continue

    count = mask.sum()
    sc.pl.embedding(
        adata[mask],
        basis="spatial",
        color=ct_col,
        size=12 if count < 2000 else 50000 / count,
        legend_loc=None,
        palette={ct: cell_type_colors.get(ct, "gray")},
        show=False,
        ax=axs[i],
        title=f"{ct} (n={count})",
    )
    axs[i].set_xlim(xlim)
    axs[i].set_ylim(ylim)
    handles.append(
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=cell_type_colors.get(ct, "gray"),
            markersize=10,
        )
    )
    labels.append(ct)

for ax in axs[len(cell_types) :]:
    fig.delaxes(ax)

fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.subplots_adjust(right=0.85)
plt.suptitle("Spatial Plots by Cell Type", fontsize=16, y=1.02)
plt.savefig(
    os.path.join(annotation_path, f"Spatial_faceted_{ct_col}.png"),
    dpi=200,
    bbox_inches="tight",
)
plt.close()

logger.info("Exporting output...")

# diagnostic
# mmc_leiden_crosstab.to_csv(
#    os.path.join(annotation_path, "mmc_leiden_crosstab_normalized.csv")
# )

# subset columns
if "cell_id" not in adata.obs.columns:
    logger.error(f"No cell_ID column. Available columns: {adata.obs.columns}")
adata.obs = adata.obs[
    [
        col
        for col in adata.obs.columns
        if any(substr in col for substr in ["cell_id", "leiden", "score", "cell_type"])
    ]
]
# replace np.nan with "None" for downstream compatibility
na_cols = [
    "cell_type_mmc_runner_up_1",
    "cell_type_mmc_runner_up_2",
    "cell_type_mmc_runner_up_1_incl_low_quality",
    "cell_type_mmc_runner_up_2_incl_low_quality",
]
for col in na_cols:
    if isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
        adata.obs[col] = adata.obs[col].cat.add_categories("None")
    adata.obs[col] = adata.obs[col].fillna("None")
adata.obs.to_csv(os.path.join(annotation_path, "adata_obs_annotated.csv"), index=False)
logger.info("Done.")
