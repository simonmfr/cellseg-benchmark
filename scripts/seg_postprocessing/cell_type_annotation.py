#!/usr/bin/env python
import datetime
import pathlib

"""Cell type annotation: MapMyCells, SCALPEL label transfer QC, cluster vote, marker revision.

1. Run MapMyCells against the ABC mouse brain atlas (or reuse an existing result)
2. SCALPEL QC (DoubleMAD per supertype, with bimodal handling): failing cells -> "Undefined"
3. Group subclasses into coarse cell types (subclasses without a group, e.g. Lymphoid -> "Undefined")
4. Leiden clustering; each cluster gets its majority label, QC-failed cells voting
   "Undefined" -> cell_type_vote
5. Marker revision with curated markers (reassign only, never "Undefined") -> cell_type_revised
6. Annotation QC summary (annotation_qc.csv), plots and adata_obs_annotated.csv
"""

import argparse
import json
import logging
import math
import os
import sys
import warnings

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from spatialdata import read_zarr

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))

import cellseg_benchmark.cell_annotation_utils as anno_utils
from cellseg_benchmark._constants import cell_type_colors
from cellseg_benchmark.dea_utils import add_ensembl_id

plt.rcParams["font.family"] = (
    "Arial" if "Arial" in [f.name for f in fm.fontManager.ttflist] else "sans-serif"
)
plt.rcParams["font.weight"] = "normal"
# see https://medium.com/@daimin0514/how-to-install-the-arial-font-in-linux-9e6ac76d3d9f

today = datetime.date.today().strftime("%Y%m%d")

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
    help="MAD_low factor (>0) for removing outlier annotations. SCALPEL uses 3.",
)
parser.add_argument(
    "--leiden_res", default=20.0, type=float, help="Leiden clustering resolution"
)
parser.add_argument(
    "--marker_min_score",
    default=1.0,
    type=float,
    help="Minimum cluster-mean marker score for a marker revision",
)
parser.add_argument(
    "--marker_delta",
    default=0.25,
    type=float,
    help="Margin over the voted label's marker score required for a marker revision",
)
args = parser.parse_args()

if args.mad_factor <= 0:
    parser.error("--mad_factor must be positive")

method_path = pathlib.Path(
    args.data_dir, "samples", args.sample_name, "results", args.seg_method
)
annotation_path = pathlib.Path(method_path, "cell_type_annotation")
os.makedirs(annotation_path, exist_ok=True)
mmc_dir = pathlib.Path(annotation_path, "mapmycells_out")
os.makedirs(mmc_dir, exist_ok=True)

logger.info("Loading data...")
adata = read_zarr(pathlib.Path(method_path, "sdata.zarr"))["table"]
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
files = [f for f in pathlib.Path(mmc_dir).glob(f"*{pattern}")]
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
        pathlib.Path(mmc_dir) / f"{today}_MapMyCells_{args.sample_name}_{args.seg_method}.json"
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

logger.info("SCALPEL QC...")
qc = anno_utils.scalpel_qc(
    allen_mmc_metadata["allen_cor_SUPT"], allen_mmc_metadata["allen_SUPT"], args.mad_factor
)
allen_mmc_metadata = allen_mmc_metadata.join(qc)
allen_mmc_metadata["allen_SUBC"] = anno_utils.group_cell_types(allen_mmc_metadata["allen_SUBC"]).fillna("Undefined")
allen_mmc_metadata["allen_SUBC_incl_low_quality"] = allen_mmc_metadata["allen_SUBC"].where(
    qc["qc_passed"], "Undefined"
)
logger.info(f"QC failed: {(~qc['qc_passed']).mean():.1%} of cells")

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
with plt.rc_context({"figure.figsize": (9, 9)}):
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
    plt.savefig(pathlib.Path(annotation_path, "UMAP_mapmycells_metrics.png"))
    plt.close()
    adata.obs.drop(columns=adata.obsm["allen_cell_type_mapping"].columns, inplace=True)

adata.obs["cell_type_mmc_is_low_quality"] = np.where(
    adata.obs["cell_type_mmc_incl_low_quality"] == "Undefined", "undefined", "mapped"
)
fig, axes = plt.subplots(1, 2, figsize=(18, 8))
for ax, basis in zip(axes, ["umap", "spatial"]):
    sc.pl.embedding(
        adata,
        basis=basis,
        color="cell_type_mmc_is_low_quality",
        size=(pt_size_umap if basis == "umap" else 150000 / adata.shape[0]),
        palette={"undefined": "blue", "mapped": "lightgrey"},
        legend_loc="right margin" if basis == "umap" else None,
        ax=ax,
        show=False,
    )
    ax.set_aspect("equal")
plt.savefig(
    pathlib.Path(annotation_path, "UMAP_and_Spatial_mapmycells_undefined.png"),
    dpi=200,
    bbox_inches="tight",
)
plt.close()
adata.obs.drop(columns="cell_type_mmc_is_low_quality", inplace=True)

leiden_col = f"leiden_res{args.leiden_res}".replace(".", "_")
if leiden_col not in adata.obs:
    sc.tl.leiden(adata, key_added=leiden_col, resolution=args.leiden_res)
adata.obs["cell_type_vote"], adata.obs["cell_type_revised"] = anno_utils.annotate_clusters(
    adata,
    cluster_col=leiden_col,
    label_col="cell_type_mmc_incl_low_quality",
    marker_csv=pathlib.Path(
        args.data_dir,
        "misc",
        "scRNAseq_ref_ABCAtlas_Yao2023Nature",
        "marker_genes_df",
        "20250416_cell_type_markers_top50.csv",
    ),
    min_score=args.marker_min_score,
    delta=args.marker_delta,
    logger=logger,
)
logger.info(f"cell_type_revised:\n{adata.obs['cell_type_revised'].value_counts()}")

qc_summary = anno_utils.annotation_qc_summary(adata)
logger.info(f"Annotation QC:\n{qc_summary.to_string()}")
pd.concat([pd.Series({"sample": args.sample_name, "seg_method": args.seg_method}), qc_summary]).to_frame().T.to_csv(
    pathlib.Path(annotation_path, "annotation_qc.csv"), index=False
)

plot_keys = [
    "cell_type_mmc_raw",
    "cell_type_mmc_incl_low_quality",
    "cell_type_vote",
    "cell_type_revised",
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

out_png = pathlib.Path(annotation_path, "UMAP_and_Spatial_annotation_results.png")
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
    pathlib.Path(annotation_path, f"Spatial_faceted_{ct_col}.png"),
    dpi=200,
    bbox_inches="tight",
)
plt.close()

logger.info("Exporting output...")

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
adata.obs.to_csv(pathlib.Path(annotation_path, "adata_obs_annotated.csv"), index=False)
logger.info("Done.")
