import math
import os
import re
import subprocess
import datetime
import io
import pathlib
import warnings

import anndata as ad
import h5py
import pandas as pd
import scanpy as sc
from anndata.io import read_elem

from .. import _constants
from ..cell_annotation_utils import group_cell_types


def read_ABCAtlas(vascular_subset=False, base_path=_constants.BASE_PATH):
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
        pathlib.Path(base_path)
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
            pathlib.Path(base_path)
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

    adata.obs["cell_type_dea"] = group_cell_types(adata.obs["subclass"]).astype("category")
    return adata


def read_adata(cohort, method=None, adata_name="adata_integrated", base_path=_constants.BASE_PATH, obs_only=False):
    """Read adata from disk.

    Args:
        cohort: cohort to use
        method: method to read adata for. If None, will choose the first method it finds in cohort
        adata_name: name of adata to read. Either adata_integrated or adata_vascular_subset
        base_path: base path to read data from.
        obs_only: read only ``.obs`` (fast, low memory), for metrics that need no expression data.

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
    if obs_only:
        with h5py.File(adata_path, "r") as f:
            return ad.AnnData(obs=read_elem(f["obs"]))
    adata = sc.read_h5ad(adata_path)
    return adata


def compute_metric_for_all_methods(
    metric_func,
    cohort,
    results_name,
    base_path=_constants.BASE_PATH,
    methods=None,
    adata_name="adata_integrated",
    overwrite=False,
    skip_3D=False,
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
        skip_3D: skip 3D methods, for metrics that only apply to 2D segmentations
        kwargs: further keyword arguments passed to metric_func

    Returns:
        The results csv as DataFrame, or None if it does not exist.
    """
    if methods is None:
        data_path = pathlib.Path(base_path) / "analysis" / cohort
        methods = os.listdir(data_path)
    if skip_3D:
        methods = [m for m in methods if not m.startswith(tuple(_constants.methods_3D))]
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
    results_path = pathlib.Path(base_path) / "metrics" / cohort / results_name
    return pd.read_csv(results_path, index_col=0) if results_path.exists() else None


def compute_metric(
    metric_func,
    cohort,
    method,
    results_name,
    adata_name="adata_integrated",
    overwrite=False,
    base_path=_constants.BASE_PATH,
    obs_only=False,
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
        obs_only: read only ``.obs`` of the adata, for metrics that need no expression data
        **kwargs: kwargs for metric_fn

    Returns:
        Nothing, saves results csv in results folder
    """
    # set up paths
    results_name = pathlib.Path(base_path) / "metrics" / cohort / results_name
    # ensure results folder exists
    results_name.parent.mkdir(parents=True, exist_ok=True)

    # check if results exist and if allowed to overwrite
    adata_path = pathlib.Path(base_path) / "analysis" / cohort / method / "adatas" / f"{adata_name}.h5ad.gz"
    if results_name.exists():
        results_df = pd.read_csv(results_name, index_col=0)
        if method in results_df["method"].unique():
            if not adata_path.exists():
                print(f"No {adata_name} for {method}, dropping stale results")
                results_df[results_df["method"] != method].to_csv(results_name)
                return
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
    adata = read_adata(cohort, method, adata_name, base_path, obs_only=obs_only)
    if adata is None:
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


def parse_slurm_mem_to_gb(x: str) -> float:
    """Helper function to parse slurm memory usage in GB."""
    if x is None:
        return math.nan

    s = str(x).strip()
    if s in ("", "Unknown", "None", "NA", "nan"):
        return math.nan
    if s == "0":
        return 0.0

    m = re.match(r"^([0-9]*\.?[0-9]+)([KMGTP]?)$", s)
    if not m:
        return math.nan

    val = float(m.group(1))
    unit = m.group(2)
    mult = {
        "": 1,
        "K": 1024,
        "M": 1024**2,
        "G": 1024**3,
        "T": 1024**4,
        "P": 1024**5,
    }[unit]
    return (val * mult) / (1024**3)


def clean_method_name(name: str) -> str:
    """Return the figure label for a results directory name."""
    if name not in _constants.method_names:
        warnings.warn(f"no figure label for method {name!r}")
    return _constants.method_names.get(name, name.replace("_", " "))


def export_job_metrics_tsv(
    ref_file_path=f"{_constants.BASE_PATH}/misc/logs/run_log.tsv",
    out_dir=f"{_constants.BASE_PATH}/misc/extracted_job_stats",
):
    """Export aggregated Slurm job metrics for all jobids in ref_file_path to:
    <out_dir>/YYYYMMDD_job_data.tsv.
    """
    legacy = pathlib.Path(ref_file_path).parent / "job_runs.tsv"
    if legacy.exists():
        raise FileNotFoundError(f"Obsolete {legacy} must be merged into run_log.tsv and deleted.")

    df = pd.read_csv(ref_file_path, sep="\t")
    df["jobid"] = df["jobid"].astype(str)

    jobids = sorted(df["jobid"].dropna().unique().tolist())
    if not jobids:
        raise ValueError("No jobids found in reference file.")

    sacct_frames = [run_sacct(jobids[i:i+200]) for i in range(0, len(jobids), 200)]

    if not sacct_frames:
        raise ValueError("No sacct data returned for the provided jobids.")

    sacct = pd.concat(sacct_frames, ignore_index=True)
    if sacct.empty:
        raise ValueError("sacct returned an empty table.")

    sacct["JobIDRaw"] = sacct["JobIDRaw"].astype(str)
    sacct["jobid_base"] = sacct["JobIDRaw"].str.split(".").str[0]
    sacct["Elapsed_s"] = pd.to_numeric(sacct["ElapsedRaw"], errors="coerce")
    sacct["AllocCPUS_n"] = pd.to_numeric(sacct["AllocCPUS"], errors="coerce")
    sacct["MaxRSS_GB"] = sacct["MaxRSS"].apply(parse_slurm_mem_to_gb)

    is_batch = sacct["JobIDRaw"].str.endswith(".batch")

    agg_all = sacct.groupby("jobid_base", as_index=False).agg(
        sacct_state=("State", "first"),
        sacct_exitcode=("ExitCode", "first"),
        elapsed_s=("Elapsed_s", "max"),
        alloccpus=("AllocCPUS_n", "max"),
    )

    agg_batch = (
        sacct[is_batch]
        .groupby("jobid_base", as_index=False)
        .agg(
            maxrss_gb_batch=("MaxRSS_GB", "max"),
        )
    )

    agg_any = sacct.groupby("jobid_base", as_index=False).agg(
        maxrss_gb_any=("MaxRSS_GB", "max"),
    )

    agg = agg_all.merge(agg_batch, on="jobid_base", how="left")
    agg = agg.merge(agg_any, on="jobid_base", how="left")
    agg["maxrss_gb"] = agg["maxrss_gb_batch"].combine_first(agg["maxrss_gb_any"])
    agg = agg.drop(columns=["maxrss_gb_batch", "maxrss_gb_any"])

    agg = agg.rename(columns={"jobid_base": "jobid"})

    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / f"{datetime.datetime.now().strftime('%Y%m%d')}_job_data.tsv"
    agg.to_csv(out_path, sep="\t", index=False)
    return out_path


def run_sacct(jobids):
    """Execute sacct to get job statistics."""
    if not jobids:
        return pd.DataFrame(
            columns=[
                "JobIDRaw",
                "State",
                "ExitCode",
                "ElapsedRaw",
                "AllocCPUS",
                "MaxRSS",
            ]
        )

    cmd = [
        "sacct",
        "--jobs",
        ",".join(map(str, jobids)),
        "--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,MaxRSS",
        "--parsable2",
        "--noheader",
    ]

    res = subprocess.run(cmd, capture_output=True, text=True, check=True)

    if not res.stdout.strip():
        return pd.DataFrame(
            columns=[
                "JobIDRaw",
                "State",
                "ExitCode",
                "ElapsedRaw",
                "AllocCPUS",
                "MaxRSS",
            ]
        )

    return pd.read_csv(
        io.StringIO(res.stdout),
        sep="|",
        header=None,
        names=["JobIDRaw", "State", "ExitCode", "ElapsedRaw", "AllocCPUS", "MaxRSS"],
        dtype=str,
    )
