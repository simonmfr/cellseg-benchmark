import os
import pathlib

import pandas as pd
import scanpy as sc
import sklearn.metrics

from . import constants


def compute_clustering_scores(
    adata, celltype_name, sample_size=1000, n_pcs=30, leiden_resolution=1, **kwargs
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
        ch = sklearn.metrics.calinski_harabasz_score(
            cur_adata.obsm["X_pca"], cur_adata.obs[celltype_name]
        )
        sh = sklearn.metrics.silhouette_score(cur_adata.obsm["X_pca"], cur_adata.obs[celltype_name])
        results.loc[sample] = {
            "calinski_harabasz_score": float(ch),
            "silhouette_score": float(sh),
        }
    # compute for all samples together
    ch = sklearn.metrics.calinski_harabasz_score(adata.obsm["X_pca"], adata.obs[celltype_name])
    sh = sklearn.metrics.silhouette_score(adata.obsm["X_pca"], adata.obs[celltype_name])
    results.loc["all"] = {
        "calinski_harabasz_score": float(ch),
        "silhouette_score": float(sh),
    }
    results = results.fillna(0)
    return results

def compute_metric_for_all_methods(
    metric_func,
    cohort,
    results_name,
    base_path=constants.BASE_PATH,
    methods=None,
    adata_name="adata_integrated",
    overwrite=False,
    pass_method=False,
    **kwargs,
):
    """Generic function to iterate over all methods to compute metric.

    Calls compute_metric on metric_func. metric_func has a mandatory argument adata, and will get passed any other kwargs.

    Args:
        metric_func: function to call for every method. Called with **kwargs.
        cohort: cohort to compute metric for
        results_name: name of the csv (and possibly subfolders) to write results to.
        base_path: base path to the data
        methods: list of methods to iterate over. If None, will iterate over all methods
        adata_name: name of adata file to read (either adata_integrated or adata_vascular_subset)
        overwrite: whether to overwrite already existing results
        pass_method: whether to pass `method` and `base_path` to metric_func (default false)
        kwargs: further keyword arguments passed to metric_func
    """
    if methods is None:
        data_path = pathlib.Path(base_path) / "analysis" / cohort
        methods = os.listdir(data_path)
    for i, method in enumerate(methods):
        print(f"{i + 1}/{len(methods)} ", end="")
        compute_metric(
            metric_func,
            cohort=cohort,
            method=method,
            results_name=results_name,
            adata_name=adata_name,
            overwrite=overwrite,
            base_path=base_path,
            pass_method=pass_method,
            **kwargs,
        )


def compute_metric(
    metric_func,
    cohort,
    method,
    results_name,
    adata_name="adata_integrated",
    overwrite=False,
    base_path=constants.BASE_PATH,
    pass_method=False,
    **kwargs,
):
    """General function to compute metric and save results as csv.

    This function does the following:
    - set up results_path in "metrics"/cohort/"results_folder"
    - check if results already exists, and manage overwriting method-specific results in the csv file (using column "method" in the csv).
    - compute metric using `metric_fn` and `**kwargs`
    - save results to csv

    Args:
        metric_func: function that computes metric / score
        cohort: cohort to compute metric for
        method: method for compute metric for
        results_name: name of the csv (and possibly subfolders) to write results to.
        adata_name: name of adata file to read (either adata_integrated or adata_vascular_subset)
        overwrite: whether to overwrite already existing results
        base_path: base path to the data
        pass_method: whether to pass method and base_path to metric_func (default false)
        **kwargs: kwargs for metric_fn

    Returns:
        Nothing, saves results csv in results folder
    """
    # set up paths
    results_name = pathlib.Path(base_path) / "metrics" / cohort / results_name
    # ensure results folder exists
    results_name.parent.mkdir(parents=True, exist_ok=True)

    # check if results exist and if allowed to overwrite
    if results_name.exists():
        results_df = pd.read_csv(results_name, index_col=0)
        if method in results_df["method"].unique():
            if not overwrite:
                print(
                    f"{metric_func.__name__} already computed for {method}. Set overwrite=True to recompute"
                )
                return
            else:
                print(f"overwriting existing results for {metric_func.__name__}")
                # remove rows with this method to overwrite with new results
                results_df = results_df[results_df["method"] != method]
    else:
        results_df = None

    print(f"Running {metric_func.__name__} for {method}")
    # read adata
    adata = read_adata(cohort, method, adata_name, base_path)
    if adata is None:
        # read error ocurred, return
        return

    # compute metric
    results = metric_func(adata, method=method, base_path=base_path, **kwargs)
    if results is None:
        # some error ocurred during metric_fn, return
        return

    # add to results_df
    results.insert(loc=0, column="method", value=method)
    if results_df is None:
        results_df = results
    else:
        results_df = pd.concat([results_df, results], ignore_index=True)
    # save results
    results_df.to_csv(results_name)


def read_adata(cohort, method=None, adata_name="adata_integrated", base_path=constants.BASE_PATH):
    """Read adata from disk.

    Args:
        cohort: cohort to use
        method: method to read adata for. If None, will choose the first method it finds in cohort
        adata_name: name of adata to read. Either adata_integrated or adata_vascular_subset
        base_path: base path to read data from.

    Returns:
        adata or None if it does not exist
    """
    if method is None:
        methods = os.listdir(pathlib.Path(base_path) / "analysis" / cohort)
        method = methods[0]
        print(f"Reading {adata_name} for method {method}")
    data_path = pathlib.Path(base_path) / "analysis" / cohort / method
    adata_path = data_path / "adatas" / f"{adata_name}.h5ad.gz"
    if not adata_path.exists():
        print(f"No adata found for cohort {cohort}, method {method}, name {adata_name}")
        return
    adata = sc.read_h5ad(adata_path)
    return adata