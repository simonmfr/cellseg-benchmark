import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse

from cellseg_benchmark import BASE_PATH
from cellseg_benchmark._constants import cell_type_colors, method_colors
from cellseg_benchmark.metrics.utils import read_ABCAtlas, read_adata

# --- Functions to compute marker genes ---


def compute_negative_markers_from_reference(
    adata, celltypes, genes, celltype_name="cell_type_dea"
):
    """Compute negative markers for single-cell reference data.

    Adapted from: https://github.com/Moldia/Xenium_benchmarking/blob/1908b7f1f7ccbb69b8bd4eea007eb3db2e0c06b5/notebooks/5_segmentation_benchmark/metrics.py#L109

    Args:
        adata: single cell anndata to compute markers from
        celltypes: list of celltypes to compute markers for
        genes: list of genes to compute markers for
        celltype_name: name of obs column to use for celltypes

    Returns:
        pd.DataFrames neg_marker_mask and ratio_celltype with shape celltypes x genes,
        containing negative markers and ratio of cells expressing each gene within a celltype
    """
    # Set threshold parameters
    min_number_cells = 10  # minimum number of cells belonging to a cluster to consider it in the analysis
    max_ratio_cells = 0.005  # maximum ratio of cells expressing a marker to call it a negative marker gene-ct pair

    # Subset adata to genes of interest (measured in spatial data)
    adata = adata[:, genes]
    # Get cell types that we find in both modalities
    shared_celltypes = list(
        set(celltypes).intersection(adata.obs[celltype_name].unique())
    )

    # Filter cell types by minimum number of cells
    celltype_count = adata.obs[celltype_name].value_counts().loc[shared_celltypes]
    ct_filter = celltype_count >= min_number_cells
    celltypes = celltype_count.loc[ct_filter].index.tolist()

    # Filter cells to eligible cell types
    adata = adata[adata.obs[celltype_name].isin(celltypes)]

    # get ratio of positive cells per cell type
    count_per_ct = sc.get.aggregate(adata, celltype_name, "count_nonzero")
    count_per_ct.obs["count"] = adata.obs[celltype_name].value_counts()
    ratio_celltype = (
        count_per_ct.layers["count_nonzero"]
        / np.array(count_per_ct.obs["count"])[:, np.newaxis]
    )

    # Get gene-cell type pairs with negative marker expression
    neg_marker_mask = np.array(ratio_celltype < max_ratio_cells)
    print(
        "Number of negative markers: ",
        {
            key: int(value)
            for key, value in zip(
                count_per_ct.obs[celltype_name], neg_marker_mask.sum(axis=1)
            )
        },
    )

    neg_marker_mask_df = pd.DataFrame(
        data=neg_marker_mask,
        index=count_per_ct.obs[celltype_name],
        columns=count_per_ct.var_names,
    )
    ratio_celltype_df = pd.DataFrame(
        data=ratio_celltype,
        index=count_per_ct.obs[celltype_name],
        columns=count_per_ct.var_names,
    )

    return neg_marker_mask_df, ratio_celltype_df


def get_negative_markers(
    cohort, vascular_subset=False, overwrite=False, base_path=BASE_PATH
):
    """Get or compute negative markers.

    Reads/writes negative markers to misc/scRNAseq_ref_ABCAtlas_Yao2023Nature/negative_markers.

    Args:
        cohort: cohort to get markers for
        vascular_subset: get markers for vascular subset
        overwrite: recompute markers and save them again
        base_path: base path to the data
    Returns:
        pd.DataFrames neg_marker_mask and ratio_celltype with shape celltypes x genes,
        containing negative markers and ratio of cells expressing each gene within a celltype
    """
    data_path = (
        Path(base_path)
        / "misc"
        / "scRNAseq_ref_ABCAtlas_Yao2023Nature"
        / "negative_markers"
    )
    data_path.mkdir(parents=True, exist_ok=True)
    marker_fname = (
        data_path
        / f"negative_markers_{cohort}_{'vascular_subset' if vascular_subset else 'all'}.csv"
    )
    ratio_fname = (
        data_path
        / f"expr_ratio_{cohort}_{'vascular_subset' if vascular_subset else 'all'}.csv"
    )
    if overwrite or not (marker_fname.exists() and ratio_fname.exists()):
        # need to recompute markers
        print("Markers not found or overwrite set to True: computing negative markers")
        adata_sc = read_ABCAtlas(vascular_subset=vascular_subset, base_path=base_path)
        adata_sp = read_adata(
            cohort,
            method=None,
            adata_name=f"adata_{'vascular_subset' if vascular_subset else 'integrated'}",
            base_path=base_path,
        )
        genes = adata_sp.var_names
        celltypes = adata_sp.obs["cell_type_revised"].unique()
        neg_marker_mask, ratio_celltype = compute_negative_markers_from_reference(
            adata_sc, celltypes, genes, celltype_name="cell_type_dea"
        )

        # save results
        neg_marker_mask.to_csv(marker_fname)
        ratio_celltype.to_csv(ratio_fname)
    else:
        neg_marker_mask = pd.read_csv(marker_fname, index_col=0)
        ratio_celltype = pd.read_csv(ratio_fname, index_col=0)
    return neg_marker_mask, ratio_celltype


# --- Functions to read previously saved marker genes ---


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


# --- Compute & plotting functions for MECR, F1 and negative marker purity


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


def compute_marker_F1_score(
    adata, celltype_name, layer="volume_log1p_norm", threshold=1
):
    """Compute F1 scores for matching marker–cell type pairs.

    Args:
        adata: AnnData object
        celltype_name: column in adata containing celltypes
        layer: layer of gene expression in adata to use for computing marker f1 score
        threshold: expression threshold

    Returns:
        pd.DataFrame with F1 scores and related metrics
    """
    marker_dict = load_marker_gene_dict(
        subset_genes=adata.var_names, subset_celltypes=adata.obs[celltype_name].unique()
    )

    # Use the specified layer or adata.X
    X = adata.layers[layer] if layer else adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)

    obs_cell_types = adata.obs[celltype_name].astype(str).values
    present_cell_types = set(obs_cell_types)
    marker_cell_types = set(marker_dict)

    unmatched = present_cell_types - marker_cell_types
    if unmatched:
        print(
            f" Note: {len(unmatched)} cell type(s) from adata not in marker_dict: {sorted(unmatched)}"
        )

    results = []

    for cell_type, markers in marker_dict.items():
        if cell_type not in present_cell_types:
            continue  # skip if cell type not in adata

        target_mask = obs_cell_types == cell_type
        other_mask = ~target_mask

        for gene in markers:
            if gene not in adata.var_names:
                continue

            gene_idx = adata.var_names.get_loc(gene)
            expr = X[:, gene_idx].toarray().ravel()
            expressed = expr > threshold

            tp = np.sum(target_mask & expressed)
            fp = np.sum(other_mask & expressed)
            fn = np.sum(target_mask & ~expressed)

            precision = tp / (tp + fp) if (tp + fp) else 0
            recall = tp / (tp + fn) if (tp + fn) else 0
            f1 = (
                2 * precision * recall / (precision + recall)
                if (precision + recall)
                else 0
            )

            results.append(
                {
                    "cell_type": cell_type,
                    "marker_gene": gene,
                    "f1_score": f1,
                    "precision": precision,
                    "recall": recall,
                    "tp": tp,
                    "fp": fp,
                    "fn": fn,
                }
            )

    return pd.DataFrame(results)


def plot_marker_F1_score(cohort, results_suffix, show=False):
    """Plot marker F1 scores.

    Plots one plot per method showing cell-type specific F1 scores, and one plot summarising scores for each method
    """
    results_file = (
        Path(BASE_PATH)
        / "metrics"
        / cohort
        / "marker_gene_metrics"
        / f"marker_f1_score_{results_suffix}.csv"
    )
    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    results_df = pd.read_csv(results_file, index_col=0)

    # plot boxplot per method show per-celltype F1 results
    for method, results in results_df.groupby("method"):
        order = results.groupby("cell_type")["f1_score"].mean().sort_values().index
        pal = {ct: cell_type_colors[ct] for ct in order}
        fig = plt.figure(figsize=(10, 5))
        sns.boxplot(
            data=results,
            x="cell_type",
            y="f1_score",
            hue="cell_type",
            palette=pal,
            legend=True,
        )
        plt.xticks(rotation=45, ha="right")
        plt.ylabel("F1 Score")
        plt.ylim(0, 1)
        plt.title(f"F1 Score Distribution per Cell Type for {method}")
        plt.tight_layout()

        if show:
            plt.show()
        fig.savefig(
            plot_path / f"marker_f1_score_{method}_{results_suffix}.png",
            bbox_inches="tight",
        )
        plt.show()

    # plot summary boxplot showing per-method results
    # compute mean f1 score by computing mean f1 per method and cell type
    # and then mean of cell-type f1 scores to get one score per method
    mean_results = (
        results_df.groupby(["method", "cell_type"]).mean("f1_score").reset_index()
    )
    order = (
        mean_results.groupby("method").mean("f1_score").sort_values("f1_score").index
    )

    pal = {method: method_colors[method] for method in order}
    fig = plt.figure(figsize=(10, 5))
    sns.boxplot(
        data=mean_results,
        x="method",
        y="f1_score",
        hue="method",
        palette=pal,
        legend=True,
    )
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("F1 Score")
    plt.ylim(0, 1)
    plt.title("F1 Score Distribution per Method")
    plt.tight_layout()

    if show:
        plt.show()
    fig.savefig(
        plot_path / f"marker_f1_score_{results_suffix}.png", bbox_inches="tight"
    )
    plt.show()


def compute_negative_marker_purity(
    adata,
    neg_marker_mask_sc,
    ratio_celltype_sc,
    celltype_name="cell_type_revised",
    **kwargs,
):
    """Negative marker purity aims to measure read leakeage between cells in spatial datasets.

    For this, we calculate the increase in positive cells assigned in spatial datasets to pairs of genes-celltyes with
    no/very low expression in scRNAseq.

    Adapted from: https://github.com/Moldia/Xenium_benchmarking/blob/1908b7f1f7ccbb69b8bd4eea007eb3db2e0c06b5/notebooks/5_segmentation_benchmark/metrics.py#L109

    Args:
    adata : AnnData
        Annotated ``AnnData`` object with counts from spatial data
    neg_marker_mask_sc: pandas dataframe containing negative markers from single cell reference (from get_negative_markers)
    ratio_celltype_sc: pandata dataframe containing ratio of cells expressing each gene within a celltype (from get_negative_markers)
    celltype_name: name of obs column to use for celltypes
    kwargs: unused, just here to catch arguments passed automatically when using this fn with `compute_metric_for_all_methods`

    Returns:
    negative marker purity:
        Increase in proportion of positive cells assigned in spatial data to pairs of genes-celltyes with no/very low expression in scRNAseq
    """
    # Set threshold parameters - same as used in get_negative_markers
    min_number_cells = 10  # minimum number of cells belonging to a cluster to consider it in the analysis

    shared_celltypes = list(neg_marker_mask_sc.index)
    # Filter cell types by minimum number of cells
    celltype_count = adata.obs[celltype_name].value_counts().loc[shared_celltypes]
    ct_filter = celltype_count >= min_number_cells
    celltypes = celltype_count.loc[ct_filter].index.tolist()

    # Filter cells to eligible cell types
    adata = adata[adata.obs[celltype_name].isin(celltypes)]

    # get ratio of positive cells per cell type
    count_per_ct = sc.get.aggregate(
        adata, celltype_name, "count_nonzero", layer="counts"
    )
    count_per_ct.obs["count"] = adata.obs[celltype_name].value_counts()
    ratio_celltype_sp = (
        count_per_ct.layers["count_nonzero"]
        / np.array(count_per_ct.obs["count"])[:, np.newaxis]
    )
    ratio_celltype_sp = pd.DataFrame(
        data=ratio_celltype_sp,
        index=count_per_ct.obs[celltype_name],
        columns=count_per_ct.var_names,
    )

    # ensure consistent celltypes
    neg_marker_mask_sc = neg_marker_mask_sc.loc[ratio_celltype_sp.index]
    ratio_celltype_sc = ratio_celltype_sc.loc[ratio_celltype_sp.index]

    # Get pos cell ratios in negative marker-cell type pairs
    lowvals_sc = np.full_like(neg_marker_mask_sc, np.nan, dtype=np.float32)
    lowvals_sc = ratio_celltype_sc.values[neg_marker_mask_sc]
    lowvals_sp = ratio_celltype_sp.values[neg_marker_mask_sc]

    # Take the mean over the normalized expressions of the genes' negative cell types
    mean_sc_low_ratio = np.nanmean(lowvals_sc)
    mean_sp_low_ratio = np.nanmean(lowvals_sp)

    # Calculate summary metric
    negative_marker_purity = 1
    if mean_sp_low_ratio > mean_sc_low_ratio:
        negative_marker_purity -= mean_sp_low_ratio - mean_sc_low_ratio

    return pd.DataFrame({"negative_marker_purity": [negative_marker_purity]})
