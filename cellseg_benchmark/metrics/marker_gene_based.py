import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from cellseg_benchmark import BASE_PATH
from cellseg_benchmark._constants import method_colors


def load_marker_gene_dict(subset_genes=None, subset_celltypes=None):
    """Load dict of marker genes, filtered by gene and celltypes.

    From Allen mouse brain scRNA-seq atlas (Yao 2023 Nature, 4M cells)
    Computed using edgeR of pseudobulks, using author-derived cell annotations (see separate script)

    Args:
        subset_genes: list of genes to subset gene dict to
        subset_celltypes: list of celltypes to subset gene dict to

    Returns:
        marker_gene_dict, dictionary with celltypes as keys and marker genes as values.
    """
    ABCAtlas_marker_df = pd.read_csv(
        os.path.join(
            BASE_PATH,
            "misc",
            "scRNAseq_ref_ABCAtlas_Yao2023Nature",
            "marker_genes_df",
            "20250211_cell_type_markers_top15_specific.csv",
        )
    )
    # turn ABCAtlas_marker_df to dict
    cell_types = ABCAtlas_marker_df.columns.tolist()
    cell_type_dict = {}
    for cell_type in cell_types:
        # Skip the index column (0) if present
        if cell_type == "0":
            continue
        # Skip if celltype not in subset_celltypes
        if (subset_celltypes is not None) and (cell_type not in subset_celltypes):
            continue
        # Get values from column, excluding the header row
        genes = ABCAtlas_marker_df[cell_type].iloc[0:].tolist()
        # Remove any NaN values
        genes = [gene for gene in genes if pd.notna(gene)]
        if subset_genes is not None:
            genes = [gene for gene in genes if gene in subset_genes]
        cell_type_dict[cell_type] = genes
    return cell_type_dict


def _MECR_score(adata, gene_pairs, layer=None):
    """Compute the Mutually Exclusive Co-expression Rate (MECR) for each gene pair in an AnnData object.

    Adapted from: https://elihei2.github.io/segger_dev/api/validation/#segger.validation.compute_MECR
    Args:
    - adata: AnnData
        Annotated data object containing gene expression data.
    - gene_pairs: List[Tuple[str, str]]
        List of tuples representing gene pairs to evaluate.

    Returns:
    - results DataFrame
    """
    results = []
    gene_expression = adata.to_df(layer=layer)
    for gene1, gene2 in gene_pairs:
        expr_gene1 = gene_expression[gene1] > 0
        expr_gene2 = gene_expression[gene2] > 0
        both_expressed = (expr_gene1 & expr_gene2).mean()
        at_least_one_expressed = (expr_gene1 | expr_gene2).mean()
        mecr = (
            both_expressed / at_least_one_expressed if at_least_one_expressed > 0 else 0
        )
        results.append({"gene1": gene1, "gene2": gene2, "MECR": mecr})
    return pd.DataFrame(results)


def compute_MECR_score(
    adata, subset_vascular_celltypes=False, layer="volume_log1p_norm"
):
    """Compute MECR Score on specific gene pairs.

    Gene pairs are loaded using load_marker_gene_dict()

    Args:
        adata: anndata to compute MECR with
        subset_vascular_celltypes: if True, subsets gene pairs to vascular celltypes.
        layer: layer of gene expression in adata to use for computing MECR score
    Returns:
        results DataFrame
    """
    subset_celltypes = None
    if subset_vascular_celltypes:
        subset_celltypes = ["ECs", "Pericytes", "SMCs", "VLMCs"]
    marker_gene_dict = load_marker_gene_dict(
        subset_genes=adata.var_names, subset_celltypes=subset_celltypes
    )
    # process marker gene dict to gene-pair list
    # get gene-pair list for MECR
    # all with all others
    gene_pairs = []
    cell_types = list(marker_gene_dict.keys())
    for i, cell_type1 in enumerate(cell_types):
        for cell_type2 in cell_types[i + 1 :]:
            for gene1 in marker_gene_dict[cell_type1]:
                for gene2 in marker_gene_dict[cell_type2]:
                    gene_pairs.append((gene1, gene2))

    results = _MECR_score(adata, gene_pairs, layer=layer)
    return results


def plot_MECR_score(cohort, results_suffix, percentile=97, show=False):
    """Plot violin plot of MECR scores."""
    results_file = (
        Path(BASE_PATH)
        / "metrics"
        / cohort
        / "marker_gene_metrics"
        / f"MECR_score_{results_suffix}.csv"
    )
    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    results_df = pd.read_csv(results_file, index_col=0)

    # Remove outliers per method and prepare for plotting
    filtered = []
    for method, df in results_df.groupby("method"):
        threshold = np.percentile(df["MECR"], percentile)
        filtered.append(df[df["MECR"] <= threshold])
    results_df_filtered = pd.concat(filtered)

    # Order datasets by median MECR value
    method_order = (
        results_df_filtered.groupby("method")["MECR"].median().sort_values().index
    )

    # Create custom palette matching the dataset order
    custom_palette = {method: method_colors[method] for method in method_order}

    fig = plt.figure(figsize=(3.5, 8), dpi=300)
    plt.grid(True, alpha=0.3, zorder=0)

    # Create violin plot with quartile lines and custom colors
    sns.violinplot(
        y="method",
        x="MECR",
        data=results_df_filtered,
        order=method_order,
        inner="quartile",
        linewidth=0.7,
        zorder=2,
        palette=custom_palette,  # Use the custom palette instead of hue
        legend=False,
    )
    plt.xlim(right=0.41)
    plt.yticks(rotation=0, va="center")
    plt.ylabel("")
    plt.xlabel("MECR Score")
    plt.tight_layout()

    if show:
        plt.show()
    fig.savefig(plot_path / f"MECR_score_{results_suffix}.png", bbox_inches="tight")
    plt.show()
