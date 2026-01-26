from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from sklearn.metrics import calinski_harabasz_score, silhouette_score

from cellseg_benchmark._constants import BASE_PATH


def compute_clustering_scores(
    cohort,
    method,
    celltype_name,
    adata_name="adata_integrated",
    overwrite=True,
    base_path=BASE_PATH,
    sample_size=10000,
    n_pcs=30,
    leiden_resolution=1,
):
    """Compute CH / SH clustering scores per sample.

    Computes PCA (and optional leiden clustering) on the entire data and then iterates over samples

    Args:
        cohort: cohort to compute metric for
        method: method for compute metric for
        celltype_name: name of celltype column in adata.obs.
            If 'leiden', uses leiden clustering with leiden_resolution to determine clusters
        adata_name: name of adata file to read (either adata_integrated or adata_vascular_subset)
        overwrite: whether to overwrite already existing results
        base_path: base path to the data
        sample_size: downsamples data to sample_size to speed up computation
        n_pcs: number of principal components to use
        leiden_resolution: resolution to use for leiden clustering
    Returns:
        Nothing, saves results csv in results folder
    """
    # set up paths
    results_path = Path(base_path) / "metrics" / cohort / "cell_type_metrics"
    results_path.mkdir(parents=True, exist_ok=True)
    data_path = Path(base_path) / "analysis" / cohort / method

    # check if results exist and if allowed to overwrite
    results_name = results_path / f"clustering_score_{adata_name}_{celltype_name}.csv"
    if results_name.exists():
        results_df = pd.read_csv(results_name, index_col=0)
        if method in results_df["method"].unique():
            if not overwrite:
                print(
                    f"Metric already computed for {method}. Set overwrite=True to recompute"
                )
                return
            else:
                # remove rows with this method to overwrite with new results
                results_df = results_df[results_df["method"] != method]
    else:
        results_df = None

    print(f"Computing clustering score {method}")
    # read adata
    adata_path = data_path / "adatas" / f"{adata_name}.h5ad.gz"
    if not adata_path.exists():
        print(f"No adata found for cohort {cohort}, method {method}, name {adata_name}")
        return
    adata = sc.read_h5ad(adata_path)
    # check if celltype name exists
    if celltype_name not in adata.obs.columns and (celltype_name != "leiden"):
        print(f"{celltype_name} not found in adata.obs")
        return

    # compute clustering score per sample
    results = _clustering_score(
        adata, celltype_name, sample_size=sample_size, n_pcs=n_pcs
    )
    # add to results_df
    results = results.reset_index(names="sample")
    results.insert(loc=0, column="method", value=method)
    if results_df is None:
        results_df = results
    else:
        results_df = pd.concat([results_df, results], ignore_index=True)
    # save results
    results_df.to_csv(results_name)


def _clustering_score(adata, celltype_name, sample_size=10000, n_pcs=30):
    results = pd.DataFrame(columns=["calinski_harabasz_score", "silhouette_score"])
    # prepare adata
    sc.pp.sample(adata, n=sample_size)
    adata.obsm["X_pca"] = sc.pp.pca(adata.X, n_comps=n_pcs)
    if celltype_name == "leiden":
        sc.pp.neighbors(adata, use_rep="X_pca")
        sc.tl.leiden(adata)
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
