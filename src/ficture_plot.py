import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns


def plot_overlap_heatmap(ficture_DEGs, marker_dict, factor_i, marker_name=""):
    """Creates a heatmap showing the number of gene overlapping between FICTURE factors and custom marker genes."""
    factors = sorted(ficture_DEGs.index.astype(int))
    overlap = pd.DataFrame(
        [
            [
                len(
                    set(marker_dict[ct]).intersection(
                        ficture_DEGs["TopGene_fc"].iloc[i].split(", ")
                    )
                )
                for ct in marker_dict
            ]
            for i in range(len(factors))
        ],
        columns=marker_dict.keys(),
        index=factors,
    ).T

    annot = overlap.mask(overlap == 0, "")
    plt.figure(figsize=(len(factors) * 0.3, len(marker_dict) * 0.3))
    sns.heatmap(overlap, annot=annot, fmt="", annot_kws={"size": 6})
    plt.xlabel(f"{factor_i} factors")
    plt.xticks(rotation=90, fontsize=8, va="top")
    plt.title(
        f"Gene Overlap\nFactors: {factor_i}\nMarker genes: {marker_name}", fontsize=10
    )
    return plt.gcf()


def plot_factors(sdata, name, factor_i, factor_no, save=None):
    """Plots the transcript count per FICTURE factor."""
    plt.figure()
    factor_count = (
        sdata[f"{name}_all_transcript_factors"]
        .groupby(f"{factor_i}_factors")
        .count()
        .compute()["x"]
    )
    ax = sns.barplot(factor_count[0:factor_no])
    ax.set(
        xlabel=f"{factor_i} factors",
        ylabel="count",
        title="Transcript count per FICTURE factor",
    )
    plt.xticks(rotation=90, fontsize=8, va="top")
    plt.figtext(
        0.11,
        0,
        f"{factor_count[factor_no] / factor_count.sum() * 100:.2f}% of transcripts not annotated",
        fontsize=8,
        alpha=0.7,
    )
    return plt.gcf()


def sort_and_score(adata, marker_genes):
    """Cell-level scoring of marker genes in adata and filtering those present in the dataset."""
    marker_genes_in_data = {}
    for ct, markers in marker_genes.items():
        markers_found = []
        for marker in markers:
            if marker in adata.var.index:
                markers_found.append(marker)
        if markers_found:
            marker_genes_in_data[ct] = markers_found

    for key in marker_genes_in_data.keys():
        sc.tl.score_genes(adata, marker_genes_in_data[key], score_name=key)

    return adata, marker_genes_in_data


def plot_clustermap_with_marker_genes(
    adata, cell_types, factor_i, marker_name="", cmap="coolwarm", center=0, z_score=True
):
    """Plots a clustered heatmap of marker gene scores."""
    adata, marker_genes_in_data = sort_and_score(adata, cell_types)
    selected_columns = [
        col for col in adata.obs.columns if col in marker_genes_in_data.keys()
    ]
    adata_sub_obs = adata.obs[selected_columns].T.copy()
    figsize = (
        max(8, len(adata_sub_obs.columns) * 0.425),
        max(8, len(adata_sub_obs) * 0.275),
    )
    g = sns.clustermap(
        adata_sub_obs, cmap=cmap, center=center, z_score=z_score, figsize=figsize
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left")
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(), rotation=90, ha="center"
    )  # Rotate x-tick labels

    plt.title(
        f"Marker Gene Scores\nFactors: {factor_i}\nMarker genes: {marker_name}",
        fontsize=10,
        loc="left",
    )
    plt.tight_layout()
    return plt.gcf()
