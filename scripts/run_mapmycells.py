"""
Run automated cell type annotation using MapMyCells (MMC) and revise output.

1) Run MMC on mouse brain reference atlas (ABCAtlas)
2) QCs MMC output, defined low-quality mappings as "Undefined", and mixed cell type identities as "Mixed", then merged labels into adata
3) To de-noise data, run Leiden clustering and assign cell types to clusters
4) Revise annotation using marker gene scores from reference atlas
5) UMAP and Spatial plots
6) Export results as csv
"""
import math, os
import sys
import warnings
import json
from os.path import join
from pathlib import Path
from subprocess import run
import pandas as pd
import numpy as np
import scanpy as sc
from spatialdata import read_zarr

import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import matplotlib.font_manager as fm
plt.rcParams["font.family"] = "Arial" if "Arial" in [f.name for f in fm.fontManager.ttflist] else "sans-serif"
plt.rcParams['font.weight'] = 'normal'
# see https://medium.com/@daimin0514/how-to-install-the-arial-font-in-linux-9e6ac76d3d9f
import seaborn as sns

sys.path.insert(1, join(Path(__file__).parent.parent.resolve(), "src"))
import cell_annotation_utils as anno_utils

from datetime import date
today = date.today().strftime("%Y%m%d")

if len(sys.argv) < 4:
    print(f"Error: Missing required arguments. Usage: python {os.path.basename(__file__)} sample_name method_name data_dir [mad_factor]")
    sys.exit(1)

if "SLURM_CPUS_PER_TASK" in os.environ:
    sc.settings.n_jobs = int(os.environ["SLURM_CPUS_PER_TASK"])
    print(f"Using {sc.settings.n_jobs} CPU cores from SLURM allocation")

warnings.filterwarnings("ignore")

sample_name = sys.argv[1]
method_name = sys.argv[2]
data_dir = sys.argv[3]

try:
    mad_factor = int(sys.argv[4]) if int(sys.argv[4]) > 0 else 3
except IndexError:
    mad_factor = 3  # Default value if not provided

path = os.path.join(data_dir, "samples", sample_name, "results", method_name)
print("Path: "+path)
annotation_path = os.path.join(path, "cell_type_annotation")
os.makedirs(annotation_path, exist_ok=True)

adata = read_zarr(os.path.join(path, "sdata.zarr"))["table"]
# Remove Blank genes
adata = adata[:, ~adata.var_names.str.startswith('Blank')]
adata.var["gene"] = adata.var.index

# add ensembl IDs
adata = anno_utils.map_gene_symbols_to_ensembl(adata)

print(adata.var.head(3))

print(f"Running MapMyCells on mouse brain atlas...")
anno_utils.run_mapmycells(
    adata,
    sample_name=sample_name,
    method_name=method_name,
    annotation_path=annotation_path,
    data_dir=data_dir,
    today=today
)

print(f"Processing MapMyCells output...")

json_path = os.path.join(annotation_path, "mapmycells_out", f"{today}_MapMyCells_{sample_name}_{method_name}.json")
with open(json_path, "rb") as src:
        json_results = json.load(src)
allen_mmc_metadata = anno_utils.process_mapmycells_output(json_results)
print(allen_mmc_metadata.head(2))

anno_utils.plot_metric_distributions(allen_mmc_metadata, 
                                     out_path=annotation_path,
                                     file_name=f"QC_raw_metric_distributions")

# Plot MAD thresholds for SUBC and CLAS (per cell type)
for allen_key, figsize in [("SUBC", (45, 7)), ("CLAS", (13, 7))]:
        anno_utils.plot_mad_thresholds(
            allen_mmc_metadata,
            out_path=annotation_path,
            name=f"QC_mad_threshold_{allen_key}",
            group_column=f"allen_{allen_key}",
            value_column=f"allen_cor_{allen_key}",
            mad_factor=mad_factor,
            figsize=figsize,
)

print(f"Mark low-quality cells as 'Undefined'...")
# based on correlation MAD as suggested by Allen Institute
taxonomy_levels = ["CLAS", "SUBC", "SUPT", "CLUS"]
prefixes = ["allen", "allen_runner_up_1", "allen_runner_up_2"]
for level in taxonomy_levels:
    for prefix in prefixes:
        anno_utils.mark_low_quality_mappings(allen_mmc_metadata, 
                                             target_column = prefix, mad_factor = mad_factor, level = level)

# re-group cell types from *_SUBC
for prefix in prefixes:
    for suffix in ["SUBC", "SUBC_incl_low_quality"]:
        allen_mmc_metadata[f"{prefix}_{suffix}"] = anno_utils.group_cell_types(allen_mmc_metadata[f"{prefix}_{suffix}"])   

print(allen_mmc_metadata["allen_SUBC"].value_counts())

print(f"Annotate 'mixed' cells based on runner-up probability...")
mixed_df = anno_utils.create_mixed_cell_types(df = allen_mmc_metadata, 
                                              diff_threshold = 0.5)

allen_mmc_metadata = allen_mmc_metadata.merge(mixed_df, 
                                              left_index=True, 
                                              right_index=True, 
                                              how='left')
print(mixed_df.head())

# add to adata.obsm
adata.obsm["allen_cell_type_mapping"] = allen_mmc_metadata.loc[adata.obs.index]

adata = anno_utils.process_adata(adata)

# Find area log-normalized layer
area_log_key = next(
        (
            k
            for k in adata.layers.keys()
            if k.endswith("_log1p_norm") and k != "library_log1p_norm"
        ),
        None,
)

# Plot distributions
fig, axs = plt.subplots(1, 4, figsize=(20, 4), gridspec_kw={"wspace": 0.4})

titles = [
        "counts",
        "library_log1p_norm",
        f"{area_log_key}",
        "zscore (volume based)",
]
data_arrays = [
        adata.layers["counts"].sum(1),
        adata.layers["library_log1p_norm"].sum(1),
        adata.layers[area_log_key].sum(1),
        adata.layers["zscore"].sum(1),
]

for ax, data, title in zip(axs, data_arrays, titles):
        sns.histplot(data, kde=False, bins=100, ax=ax, zorder=2, alpha=1)
        ax.set_title(title, fontsize=14)
        ax.grid(True, alpha=0.2)
        if ax.get_legend():
            ax.get_legend().remove()
        for patch in ax.patches:
            patch.set_edgecolor("none")

plt.tight_layout()
plt.savefig(os.path.join(annotation_path, "QC_compare_normalizations.png"), dpi=200)
plt.close()

pt_size_umap = 220000/adata.shape[0]

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
        vmin=0, vmax=1,
        cmap="cividis",
        show=False,
    )
    plt.tight_layout()
    plt.gca().set_aspect(1)
    plt.savefig(os.path.join(annotation_path, "UMAP_mapmycells_metrics.png"))
    plt.close()
    adata.obs.drop(columns=adata.obsm["allen_cell_type_mapping"].columns, inplace=True)

# plot mixed and low-quality cells on umap and spatial plot
adata.obs["cell_type_mmc_is_low_quality"] = adata.obs["cell_type_mmc_incl_low_quality"].apply(
        lambda x: "undefined" if x == "Undefined" else "mapped"
    )
palette_dict = {
    "cell_type_mmc_is_mixed": {"mixed": "red", "unique": "lightgrey"},
    "cell_type_mmc_is_low_quality": {"undefined": "blue", "mapped": "lightgrey"},
}

fig, axes = plt.subplots(2, 2, figsize=(18, 16), gridspec_kw={'wspace': 0.001, 'hspace': 0.15})

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
        show=False
    )
    axes[0, i].set_aspect('equal')

    # Spatial
    sc.pl.embedding(
        adata,
        basis='spatial',
        color=key,
        size=150000 / adata.shape[0],
        legend_loc=None,
        palette=palette,
        ax=axes[1, i],
        show=False
    )
    axes[1, i].set_aspect('equal')

output_path = os.path.join(annotation_path, "UMAP_and_SPATIAL_mapmycells_mixed_and_undefined.png")
plt.savefig(output_path, dpi=200, bbox_inches="tight")
plt.close()

print(adata.obs['cell_type_mmc_incl_mixed'].value_counts())
print(adata.obs['cell_type_mmc_is_mixed'].value_counts())
print(adata.obs['cell_type_mmc_mixed_names'].value_counts())
print(adata.obs['cell_type_mmc_is_low_quality'].value_counts())
print(adata.obs['cell_type_mmc_incl_low_quality'].value_counts())

# Revise annotations
cell_type_colors = {
    "ECs": "#FF6464",
    "Pericytes": "#F6EC2A",
    "SMCs": "#29FBA7",
    "VLMCs": "#85B0F9",
    "ABCs": "#AEC9F5",
    "Ependymal": "#FDC000",
    "Tanycytes": "#FFE180",
    "Choroid-Plexus": "#BF9800",
    "Astrocytes": "#FE9A30",
    "Astroependymal": "#FE9A30",
    "Bergmann": "#FFD1A3",
    "Oligodendrocytes": "#4564FF",
    "OPCs": "#0095FF",
    "Microglia": "#00C088",
    "BAMs": "#00F3D8",
    "Immune-Other": "#98DF8A",
    "Neurons-Gaba": "#B449F8",
    "Neurons-Glut": "#CEB3FF",
    "Neurons-Glyc-Gaba": "#DCAEFF",
    "Neurons-Dopa": "#FCA0FF",
    #"Neurons-Immature": "#FF50E5",
    "Neurons-Granule-Immature": "#FF50E5",
    "Neurons-Other": "#FCA0FF",
    "OECs": "#9EDAE5",
    "Unknown": "#D9D9D9", # = not found in mmc dict, see process_mapmycells_output()
    "Undefined": "#D9D9D9", # = below QC threshold
    "Mixed": "#D9D9D9"
}

# Dictionary mapping cell type identifiers from MapMyCells to their corresponding Score Cell Type from sc.tl.score_genes
# Format: "MapMyCells Cell Type":"Score Cell Type"
mmc_to_score_dict = {
    "Choroid-Plexus":'Choroid-Plexus',
    "BAMs":'BAMs',
    "Ependymal":'Ependymal',
    "SMCs":'SMCs',
    "Microglia":'Microglia',
    "Pericytes":'Pericytes',
    "VLMCs":'VLMCs',
    "Immune-Other":'Immune-Other',
    "ECs":'ECs',
    "ABCs":'ABCs',
    "OPCs":'OPCs',
    "Oligodendrocytes":'Oligodendrocytes',
    "OECs":'OECs',
    "Neurons-Dopa":'Neurons-Dopa',
    "Neurons-Granule-Immature":'Neurons-Granule-Immature',
    "Astroependymal":'Astroependymal',
    "Astrocytes":'Astrocytes',
    "Neurons-Other":'Neurons-Other',
    "Tanycytes":'Tanycytes',
    "Bergmann":'Bergmann',
    "Neurons-Dopa":'Neurons-Dopa-Gaba',
    "Neurons-Gaba":'Neurons-Gaba',
    "Neurons-Glyc-Gaba":'Neurons-Glyc-Gaba',
    "Neurons-Glut":'Neurons-Glut'
}

leiden_res = 10.0
leiden_col = f"leiden_res{leiden_res}".replace(".", "_")
adata, annotation_results, normalized_percentage = anno_utils.revise_annotations(
    adata,
    leiden_res=leiden_res,
    leiden_col=leiden_col, 
    cell_type_colors=cell_type_colors, 
    mmc_to_score_dict=mmc_to_score_dict,
    score_threshold=0.5,
    top_n_genes=50,
    ABCAtlas_marker_df_path=os.path.join(
                data_dir,
                "misc",
                "scRNAseq_ref_ABCAtlas_Yao2023Nature",
                "marker_genes_df",
                "20250416_cell_type_markers_top50.csv")
)

# plot main annotations as umap and spatial plot
plot_keys = ["cell_type_mmc_raw", "cell_type_mmc_raw_clusters", "cell_type_mmc_raw_revised",
             "cell_type_mmc_incl_mixed", "cell_type_mmc_incl_mixed_clusters", "cell_type_mmc_incl_mixed_revised",
             "cell_type_mmc_incl_low_quality", "cell_type_mmc_incl_low_quality_clusters", "cell_type_mmc_incl_low_quality_revised", 
             leiden_col
            ]

# colors for large number of leiden clusters
n_clusters = len(np.unique(adata.obs[leiden_col]))
base_colors = (
    sns.color_palette("tab20", 20) + 
    sns.color_palette("tab20b", 20) + 
    sns.color_palette("tab20c", 20) + 
    sns.color_palette("husl", 30) +
    sns.color_palette("Set1", 9) +
    sns.color_palette("Set2", 8) +
    sns.color_palette("Set3", 12)
)
long_palette = (base_colors * (n_clusters // len(base_colors) + 1))[:n_clusters]

# ensure consistent category order
for key in plot_keys:
    if key != leiden_col:
        adata.obs[key] = pd.Categorical(
            adata.obs[key], categories=list(cell_type_colors.keys())
        ).remove_unused_categories()

fig, axes = plt.subplots(len(plot_keys), 3, figsize=(26, 9*len(plot_keys)), gridspec_kw={'wspace': 0, 'hspace': 0.1})

for i, key in enumerate(plot_keys):

    palette = long_palette if key == leiden_col else cell_type_colors
    
    sc.pl.umap(
        adata,
        color=key,
        size=pt_size_umap,
        legend_fontoutline=2,
        legend_fontweight="normal",
        legend_fontsize=7,
        legend_loc=None,
        palette=palette,
        ax=axes[i, 0],
        show=False
    )
    axes[i, 0].set_aspect('equal')

    sc.pl.umap(
        adata,
        color=key,
        legend_loc = "on data",
        palette=palette,
        size=pt_size_umap,
        legend_fontweight='normal',
        legend_fontsize=12,
        legend_fontoutline=1.5,
        ax=axes[i, 1],
        show=False
    )
    axes[i, 1].set_aspect('equal')

    # Spatial
    sc.pl.embedding(
        adata,
        basis='spatial',
        color=key,
        size=150000 / adata.shape[0],
        palette=palette,
        legend_loc=None if key == leiden_col else 'right margin',
        ax=axes[i, 2],
        show=False
    )
    axes[i, 2].set_aspect('equal')

output_path = os.path.join(annotation_path, "UMAP_and_SPATIAL_mapmycells_results.png")
plt.savefig(output_path, dpi=200, bbox_inches="tight")
plt.show()

# plot facet plot per cell type of main annotation
ct_col = "cell_type_mmc_raw_revised"
spatial = adata.obsm['spatial']
cell_types = adata.obs[ct_col].cat.categories
n_types = len(cell_types)
n_cols = math.ceil(math.sqrt(n_types))
n_rows = math.ceil(n_types / n_cols)

fig, axs = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
axs = axs.flatten() if n_types > 1 else [axs]

# Spatial limits with margin
x, y = spatial[:, 0], spatial[:, 1]
margin_x, margin_y = 0.05 * (x.max() - x.min()), 0.05 * (y.max() - y.min())
xlim, ylim = (x.min() - margin_x, x.max() + margin_x), (y.min() - margin_y, y.max() + margin_y)

handles, labels = [], []
for i, ct in enumerate(cell_types):
    if i >= len(axs): break
    mask = adata.obs[ct_col] == ct
    if not mask.any(): continue

    count = mask.sum()
    sc.pl.embedding(
        adata[mask], basis='spatial', color=ct_col,
        size=12 if count < 2000 else 50000 / count,
        legend_loc=None, palette={ct: cell_type_colors.get(ct, 'gray')},
        show=False, ax=axs[i], title=f"{ct} (n={count})"
    )
    axs[i].set_xlim(xlim)
    axs[i].set_ylim(ylim)
    handles.append(plt.Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=cell_type_colors.get(ct, 'gray'), markersize=10))
    labels.append(ct)

# Remove unused axes
for ax in axs[len(cell_types):]: fig.delaxes(ax)

fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.subplots_adjust(right=0.85)
plt.suptitle("Spatial Plots by Cell Type", fontsize=16, y=1.02)
plt.savefig(os.path.join(annotation_path, f"SPATIAL_subplots_{ct_col}.png"), dpi=200, bbox_inches="tight")
plt.show()

# export data
normalized_percentage.to_csv(
    os.path.join(annotation_path, "mapmycells_leiden_crosstab_normalized.csv")
)
# subset columns
adata.obs = adata.obs[[col for col in adata.obs.columns if any(substr in col for substr in ['cell_id', 'leiden', 'score', 'mmc'])]]
# replace np.nan with "None" for downstream compatibility
na_cols = ['cell_type_mmc_runner_up_1', 'cell_type_mmc_runner_up_2',
           'cell_type_mmc_runner_up_1_incl_low_quality', 'cell_type_mmc_runner_up_2_incl_low_quality']
for col in na_cols:
    if pd.api.types.is_categorical_dtype(adata.obs[col]):
        adata.obs[col] = adata.obs[col].cat.add_categories("None")
    adata.obs[col] = adata.obs[col].fillna("None")
adata.obs.to_csv(os.path.join(annotation_path, "adata_obs_annotated.csv"), index=False)
