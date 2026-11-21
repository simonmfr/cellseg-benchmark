from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import sklearn.metrics

from .. import _constants
from .utils import clean_method_name


def compute_clustering_scores(
    adata, celltype_name, sample_size=1000, n_pcs=30, leiden_resolution=1, **kwargs
):
    """Compute CH / SH clustering scores per sample.

    Computes PCA (and optional leiden clustering) on the entire data and then iterates over samples

    Args:
        adata: anndata to compute score with
        celltype_name: name of celltype column in adata.obs.
            If 'leiden', uses leiden clustering with leiden_resolution to determine clusters
        sample_size: cells per sample (random subsample) to speed up computation
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
    adata, celltype_name, sample_size=1000, n_pcs=30, leiden_resolution=1
):
    results = pd.DataFrame(columns=["calinski_harabasz_score", "silhouette_score"])
    keep = adata.obs.groupby("sample", observed=True, group_keys=False).apply(
        lambda d: d.sample(min(len(d), sample_size), random_state=0)
    ).index
    adata = adata[keep].copy()
    adata.obsm["X_pca"] = sc.pp.pca(adata.X, n_comps=n_pcs)
    if celltype_name == "leiden":
        sc.pp.neighbors(adata, use_rep="X_pca")
        sc.tl.leiden(adata, resolution=leiden_resolution)
        print(
            f"Computed leiden clustering with {len(adata.obs['leiden'].unique())} clusters"
        )
    for sample in adata.obs["sample"].unique():
        cur_adata = adata[adata.obs["sample"] == sample]
        if not 1 < cur_adata.obs[celltype_name].nunique() < cur_adata.n_obs:
            continue
        ch = sklearn.metrics.calinski_harabasz_score(
            cur_adata.obsm["X_pca"], cur_adata.obs[celltype_name]
        )
        sh = sklearn.metrics.silhouette_score(cur_adata.obsm["X_pca"], cur_adata.obs[celltype_name])
        results.loc[sample] = {
            "calinski_harabasz_score": float(ch),
            "silhouette_score": float(sh),
        }
    # compute for all samples together
    if not 1 < adata.obs[celltype_name].nunique() < adata.n_obs:
        return results
    ch = sklearn.metrics.calinski_harabasz_score(adata.obsm["X_pca"], adata.obs[celltype_name])
    sh = sklearn.metrics.silhouette_score(adata.obsm["X_pca"], adata.obs[celltype_name])
    results.loc["all"] = {
        "calinski_harabasz_score": float(ch),
        "silhouette_score": float(sh),
    }
    results = results.fillna(0)
    return results


def plot_clustering_scores(cohort, results_suffix, show=False):
    """Box plots of silhouette and Calinski-Harabasz scores per method (samples, ordered by mean silhouette) and their correlation; "merged" averages cohorts equally."""
    results_file = Path(_constants.BASE_PATH) / "metrics" / cohort / "cell_type_metrics" / f"clustering_score_{results_suffix}.csv"
    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)
    cohorts = ["aging", "ABCAtlas", "VizgenMouseBrain"] if cohort == "merged" else [cohort]
    scores = pd.concat(
        [pd.read_csv(str(results_file).replace(f"/{cohort}/", f"/{c}/"), index_col=0).assign(cohort=c) for c in cohorts],
        ignore_index=True,
    ).query("sample != 'all'")
    scores["method"] = scores["method"].map(clean_method_name)
    cols = {"silhouette_score": "Silhouette", "calinski_harabasz_score": "Calinski-Harabasz"}
    means = scores.groupby(["method", "cohort"])[list(cols)].mean().groupby("method").mean()
    order = means.sort_values("silhouette_score", ascending=False).index

    fig, axs = plt.subplots(1, 3, figsize=(9, 0.16 * len(order) + 0.8), gridspec_kw={"width_ratios": [1, 1, 1.3]})
    for ax, (col, label) in zip(axs, cols.items()):
        sns.boxplot(scores, x=col, y="method", order=order, ax=ax, color="0.85", linewidth=0.6, fliersize=1.5)
        ax.set(xlabel=label, ylabel="")
    axs[1].set_yticklabels([])
    axs[2].scatter(means.calinski_harabasz_score, means.silhouette_score, s=10, color="0.3")
    rho = means.corr("spearman").iloc[0, 1]
    axs[2].set(xlabel="Calinski-Harabasz", ylabel="Silhouette", title=f"Method means, Spearman ρ = {rho:.2f}")
    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(plot_path / f"cluster_scores_{results_suffix}.png", dpi=300)
    plt.show() if show else plt.close(fig)
