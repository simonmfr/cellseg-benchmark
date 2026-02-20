import os
from pathlib import Path

import pandas as pd
import scanpy as sc

from cellseg_benchmark import BASE_PATH


def read_ABCAtlas(vascular_subset=False, base_path=BASE_PATH):
    """Convenience function to read ABCAtlas adata.

    Args:
        vascular_subset: read vascular subset of ABCAtlas
        base_path: base path to read data from
    """
    if vascular_subset:
        adata_name = "20250205_merged_v3_vascular_subset.h5ad.gz"
    else:
        adata_name = "20250129_merged_v3.h5ad.gz"

    fname = (
        Path(base_path)
        / "misc"
        / "scRNAseq_ref_ABCAtlas_Yao2023Nature"
        / "anndata-objects"
        / adata_name
    )

    # read adata
    adata = sc.read_h5ad(fname)

    if not vascular_subset:
        # read cell metadata
        cell_meta_fname = (
            Path(base_path)
            / "misc"
            / "scRNAseq_ref_ABCAtlas_Yao2023Nature"
            / "20250411_meta_data_formatted.csv.gz"
        )
        cell_meta = pd.read_csv(
            cell_meta_fname,
            index_col=0,
            dtype=str,
        )
        cell_meta = cell_meta.astype("category")

        # join with obs
        adata.obs = adata.obs.join(cell_meta, how="left")

    return adata


def read_adata(cohort, method=None, adata_name="adata_integrated", base_path=BASE_PATH):
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
        methods = os.listdir(Path(base_path) / "analysis" / cohort)
        method = methods[0]
        print(f"Reading {adata_name} for method {method}")
    data_path = Path(base_path) / "analysis" / cohort / method
    adata_path = data_path / "adatas" / f"{adata_name}.h5ad.gz"
    if not adata_path.exists():
        print(f"No adata found for cohort {cohort}, method {method}, name {adata_name}")
        return
    adata = sc.read_h5ad(adata_path)
    return adata


def compute_metric_for_all_methods(
    metric_func,
    cohort,
    results_name,
    base_path=BASE_PATH,
    methods=None,
    adata_name="adata_integrated",
    overwrite=False,
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
        kwargs: further keyword arguments passed to metric_func
    """
    if methods is None:
        data_path = Path(base_path) / "analysis" / cohort
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
            **kwargs,
        )


def compute_metric(
    metric_func,
    cohort,
    method,
    results_name,
    adata_name="adata_integrated",
    overwrite=False,
    base_path=BASE_PATH,
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
        **kwargs: kwargs for metric_fn

    Returns:
        Nothing, saves results csv in results folder
    """
    # set up paths
    results_name = Path(base_path) / "metrics" / cohort / results_name
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
    results = metric_func(adata, **kwargs)
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
