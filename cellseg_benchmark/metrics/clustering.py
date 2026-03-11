from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from sklearn.metrics import calinski_harabasz_score, silhouette_score

from cellseg_benchmark import BASE_PATH


def compute_clustering_scores(
    adata, celltype_name, sample_size=1000, n_pcs=30, leiden_resolution=1
):
    """Compute CH / SH clustering scores per sample.

    Computes PCA (and optional leiden clustering) on the entire data and then iterates over samples

    Args:
        adata: anndata to compute score with
        celltype_name: name of celltype column in adata.obs.
            If 'leiden', uses leiden clustering with leiden_resolution to determine clusters
        sample_size: downsamples data to sample_size to speed up computation
        n_pcs: number of principal components to use
        leiden_resolution: resolution to use for leiden clustering
    Returns:
        results DataFrame or None
    """
    # check if celltype name exists
    if celltype_name not in adata.obs.columns and (celltype_name != "leiden"):
        print(f"{celltype_name} not found in adata.obs")
        return None
    # compute clustering score per sample
    results = _clustering_score(
        adata,
        celltype_name,
        sample_size=sample_size,
        n_pcs=n_pcs,
        leiden_resolution=leiden_resolution,
    )
    return results.reset_index(names="sample")


def _clustering_score(
    adata, celltype_name, sample_size=10000, n_pcs=30, leiden_resolution=1
):
    results = pd.DataFrame(columns=["calinski_harabasz_score", "silhouette_score"])
    # prepare adata
    sc.pp.sample(adata, n=sample_size)
    adata.obsm["X_pca"] = sc.pp.pca(adata.X, n_comps=n_pcs)
    if celltype_name == "leiden":
        sc.pp.neighbors(adata, use_rep="X_pca")
        sc.tl.leiden(adata, resolution=leiden_resolution)
        print(
            f"Computed leiden clustering with {len(adata.obs['leiden'].unique())} clusters"
        )
    for sample in adata.obs["sample"].unique():
        cur_adata = adata[adata.obs["sample"] == sample]
        ch = calinski_harabasz_score(
            cur_adata.obsm["X_pca"], cur_adata.obs[celltype_name]
        )
        sh = silhouette_score(cur_adata.obsm["X_pca"], cur_adata.obs[celltype_name])
        results.loc[sample] = {
            "calinski_harabasz_score": float(ch),
            "silhouette_score": float(sh),
        }
    # compute for all samples together
    ch = calinski_harabasz_score(adata.obsm["X_pca"], adata.obs[celltype_name])
    sh = silhouette_score(adata.obsm["X_pca"], adata.obs[celltype_name])
    results.loc["all"] = {
        "calinski_harabasz_score": float(ch),
        "silhouette_score": float(sh),
    }
    results = results.fillna(0)
    return results


def plot_clustering_scores(cohort, results_suffix, show=False):
    """Plot box plot of clustering scores and dot plot comparing SH and CH scores."""
    results_file = (
        Path(BASE_PATH)
        / "metrics"
        / cohort
        / "cell_type_metrics"
        / f"clustering_score_{results_suffix}.csv"
    )
    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    scores_df = pd.read_csv(results_file, index_col=0)
    # remove sample all from scores
    scores_df = scores_df[scores_df["sample"] != "all"]

    # mean scores computation
    mean_scores = scores_df.groupby("method").mean(
        ["calinski_harabasz_score", "silhouette_score"]
    )

    # Set theme
    plt.style.use("seaborn-v0_8-whitegrid")
    # --- Figure 1: boxplot of CH/SH scores ---
    palette = sns.color_palette("viridis_r", n_colors=len(mean_scores))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

    # box plot of CH scores
    ch_order = list(
        mean_scores.sort_values("calinski_harabasz_score", ascending=False).index
    )
    _ = sns.boxplot(
        scores_df,
        x="method",
        y="calinski_harabasz_score",
        order=ch_order,
        ax=ax1,
        palette=palette,
    )
    ax1.set_title("Calinski-Harabasz Score", fontsize=12)
    ax1.set_ylabel("Score (higher is better)", fontsize=12)

    # box plot of CH scores
    sh_order = list(mean_scores.sort_values("silhouette_score", ascending=False).index)
    _ = sns.boxplot(
        scores_df,
        x="method",
        y="silhouette_score",
        order=sh_order,
        ax=ax2,
        palette=palette,
    )
    ax2.set_title("Silhouette Score", fontsize=12)
    ax2.set_ylabel("Score (closer to 1 is better)", fontsize=12)

    for ax, order in [(ax1, ch_order), (ax2, sh_order)]:
        # Fix tick labels alignment - set x-ticks explicitly
        ax.set_xticks(np.arange(len(order)))
        ax.set_xticklabels(order, rotation=45, ha="right", fontweight="medium")
        ax.set_xlabel("")
        # Add subtle border
        for spine in ax.spines.values():
            spine.set_edgecolor("#dddddd")
            spine.set_linewidth(0.8)

    plt.tight_layout()
    if show:
        plt.show()
    fig.savefig(plot_path / f"cluster_scores_{results_suffix}.png")
    plt.show()

    # --- Figure 2: scatter plot comparing CH/SH scores ---
    fig = plt.figure(figsize=(7, 5))

    # Create norm for colormap
    # norm = plt.Normalize(
    #    mean_scores["calinski_harabasz_score"].min(),
    #    mean_scores["calinski_harabasz_score"].max(),
    # )

    # Scatter plot with viridis colormap based on CH score
    scatter = plt.scatter(
        mean_scores["calinski_harabasz_score"],
        mean_scores["silhouette_score"],
        c=mean_scores["calinski_harabasz_score"],
        cmap="viridis",
        s=80,
        alpha=0.8,
        edgecolor="white",
        linewidth=0.5,
    )

    # Add colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label("Calinski-Harabasz Score", fontsize=10)

    # Add labels for each point
    for i, txt in enumerate(mean_scores.index):
        plt.annotate(
            txt,
            (
                mean_scores["calinski_harabasz_score"].iloc[i],
                mean_scores["silhouette_score"].iloc[i],
            ),
            fontsize=9,
            ha="center",
            va="bottom",
            xytext=(0, 5),
            textcoords="offset points",
        )

    plt.xlabel("Calinski-Harabasz Score", fontsize=12)
    plt.ylabel("Silhouette Score (closer to 1 is better)", fontsize=12)
    plt.grid(True, alpha=0.3)

    # Add a box around the plot with light gray color
    plt.gca().spines["top"].set_visible(True)
    plt.gca().spines["right"].set_visible(True)
    for spine in plt.gca().spines.values():
        spine.set_edgecolor("#dddddd")
        spine.set_linewidth(0.8)

    plt.tight_layout()
    if show:
        plt.show()
    fig.savefig(plot_path / f"cluster_scores_dotplot_{results_suffix}.png")
    plt.show()
