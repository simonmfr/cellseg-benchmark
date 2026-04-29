import math
import os
import re
import subprocess
import datetime
import io
import pathlib
from typing import Union, Optional

import pandas as pd
import scanpy as sc

from .. import _constants


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

    return adata


def read_adata(cohort, method=None, adata_name="adata_integrated", base_path=_constants.BASE_PATH):
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


def compute_metric_for_all_methods(
    metric_func,
    cohort,
    results_name,
    base_path=_constants.BASE_PATH,
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
    base_path=_constants.BASE_PATH,
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


def normalize_jobname(jobname: str) -> str:
    """Helper function to fix typo in jobname of vpt."""
    if not isinstance(jobname, str):
        return ""
    # fix historic typo: vtp -> vpt
    return re.sub(r"^vtp", "vpt", jobname)


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


def method_with_flavor_from_row(jobname: str, key: str) -> str:
    """Create a canonical method string."""
    j = normalize_jobname(jobname)
    k = str(key)

    # --- Baysor with optional qualifier before key
    m = re.match(
        rf"^Baysor_(?:(?P<qualifier>.+?)_)?{re.escape(k)}_CP(?P<cp>\d+)_(?P<stain>[^_]+)_(?P<conf>.+)$",
        j,
    )
    if m:
        qualifier = m.group("qualifier")
        cp = m.group("cp")
        stain = m.group("stain")
        conf = m.group("conf")
        prefix = "Baysor"
        if qualifier:
            prefix += f"_{qualifier}"
        if stain == "nuclei":
            return f"{prefix}_2D_Cellpose_{cp}_nuclei_model_{conf}"
        return f"{prefix}_2D_Cellpose_{cp}_DAPI_{stain}_{conf}"

    # --- vpt2D / vpt3D: vpt2D_<key>_<stain...> -> vpt2D_<stain...>
    for dim in ("2D", "3D"):
        prefix = f"vpt{dim}_{k}_"
        if j.startswith(prefix):
            rest = j[len(prefix) :]
            return f"vpt_{dim}_DAPI_{rest}"
        if j == f"vpt{dim}_{k}":
            return f"vpt_{dim}"

    # --- Proseg_3D vpt: Proseg_3D_<key>_vpt2D_<flavor>_vxl_<voxel>
    prefix_3d = f"Proseg_3D_{k}_"
    if j.startswith(prefix_3d) and "_vpt" in j:
        rest = j[len(prefix_3d):]
        rest = re.sub(r"_vxl_.+$", "", rest)  # ignore voxel size
        m = re.match(r"^vpt(?P<dim>2D|3D)_(?P<flavor>.+)$", rest)
        if m:
            return f"Proseg_3D_vpt_{m.group('dim')}_DAPI_{m.group('flavor')}"
        return f"Proseg_3D_{rest}"

    # --- Proseg_3D CP: Proseg_3D_<key>_CP2_PolyT_vxl_7
    if j.startswith(prefix_3d) and "_CP" in j:
        rest = j[len(prefix_3d):]
        m = re.match(r"^CP(?P<cp>\d+)_(?P<stain>[^_]+)(?:_vxl_.+)?$", rest)
        if m:
            cp = m.group("cp")
            stain = m.group("stain")
            if stain == "nuclei":
                return f"Proseg_3D_Cellpose_{cp}_nuclei_model"
            return f"Proseg_3D_Cellpose_{cp}_DAPI_{stain}"
        return f"Proseg_3D_{re.sub(r'_vxl_.+$', '', rest)}"

    # --- Proseg (2D) CP: Proseg_<key>_CP2_PolyT_vxl_7
    prefix_2d = f"Proseg_{k}_"
    if j.startswith(prefix_2d) and "_CP" in j:
        rest = j[len(prefix_2d):]
        m = re.match(r"^CP(?P<cp>\d+)_(?P<stain>[^_]+)(?:_vxl_.+)?$", rest)
        if m:
            cp = m.group("cp")
            stain = m.group("stain")
            if stain == "nuclei":
                return f"Proseg_2D_Cellpose_{cp}_nuclei_model"
            return f"Proseg_2D_Cellpose_{cp}_DAPI_{stain}"
        return f"Proseg_{re.sub(r'_vxl_.+$', '', rest)}"

    # --- Cellpose jobs: CP1_<key>_<stain> / CP2_<key>_<stain>
    prefix = f"CP1_{k}_"
    if j.startswith(prefix):
        stain = j[len(prefix):]
        if stain == "nuclei":
            return "Cellpose_1_nuclei_model"
        if stain == "Merlin":
            return "Cellpose_1_Merlin"
        return f"Cellpose_1_DAPI_{stain}"

    prefix = f"CP2_{k}_"
    if j.startswith(prefix):
        stain = j[len(prefix):]
        if stain == "nuclei":
            return "Cellpose_2_nuclei_model"
        return f"Cellpose_2_DAPI_{stain}"

    # --- Simple methods
    if j == f"voronoi_{k}":
        return "Negative_Control_Voronoi"
    if j == f"nuclei_{k}":
        return "Negative_Control_Nuclei"
    if j == f"merscope_{k}":
        return "MERSCOPE"
    if j == f"transcript_tif_{k}":
        return "Transcript_TIF"
    m = re.match(rf"^rastered(?P<width>\d+)_{re.escape(k)}$", j)
    if m:
        return f"Negative_Control_Rastered_{m.group('width')}"

    return j


def find_latest_job_data_tsv(metrics_dir):
    """Return the newest exported job-data TSV in metrics_dir.
    Uses file modification time and only considers '*_job_data.tsv'.
    """
    metrics_dir = pathlib.Path(metrics_dir)
    candidates = [p for p in metrics_dir.glob("*_job_data.tsv") if p.is_file()]
    if not candidates:
        raise FileNotFoundError(f"No '*_job_data.tsv' files found in: {metrics_dir}")
    return max(candidates, key=lambda p: p.stat().st_mtime)


def export_job_metrics_tsv(
    ref_file_path="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/job_runs.tsv",
    out_dir="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/extracted_job_stats",
):
    """Export aggregated Slurm job metrics for all jobids in ref_file_path to:
    <out_dir>/YYYYMMDD_job_data.tsv.
    """
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


def prepare_label_maps(factor_to_celltype, true_cluster):
    """Prepare label maps."""
    correct_celltypes = {
        factor_to_celltype[i]: true_cluster[factor_to_celltype[i]]
        for i in factor_to_celltype.keys()
    }

    factor_map = factor_to_celltype.copy()
    factor_map["-1"] = "unassigned"

    return factor_map, correct_celltypes

def show_all_method_names(
        ref_file_path: Union[str, pathlib.Path],
        only_successful: bool=False,
        metrics_dir: Optional[Union[str, pathlib.Path]]=None,
        return_counts: bool=False
):
    """
    Show all possible canonical method names derived from the ref file.

    Args:
        ref_file_path (str | Path): Path to the appended job_runs.tsv file.
        only_successful (bool): If True, only include methods with successful runs according to the
            newest *_job_data.tsv file in metrics_dir. Defaults to False.
        metrics_dir (str | Path, optional): Required if only_successful=True. Path to the metrics directory of sacct.
        return_counts (bool): If True, return a DataFrame with counts instead of a plain list. Defaults to False.

    Returns:
        list[str] or pd.DataFrame: List of all method names in the job_data.tsv file.

    Note:
        If only_successful=True requires prior run of export_job_metrics_tsv() with access to sacct-command.
    """

    df = pd.read_csv(ref_file_path, sep="\t")
    df["_ref_order"] = range(len(df))

    df["jobid"] = df["jobid"].astype(str)
    df["jobname"] = df["jobname"].astype(str)
    df["sample"] = df["key"].astype(str)
    df["jobname_norm"] = df["jobname"].apply(normalize_jobname)

    df["method_with_flavor"] = df.apply(
        lambda r: method_with_flavor_from_row(r["jobname"], r["sample"]),
        axis=1,
    )

    if only_successful:
        if metrics_dir is None:
            raise ValueError("metrics_dir must be provided when only_successful=True.")

        latest_metrics_file = find_latest_job_data_tsv(metrics_dir)
        sacct = pd.read_csv(latest_metrics_file, sep="\t")
        sacct["jobid"] = sacct["jobid"].astype(str)

        required_cols = {"jobid", "sacct_state", "sacct_exitcode"}
        missing = required_cols.difference(sacct.columns)
        if missing:
            raise ValueError(
                f"Metrics file {latest_metrics_file} is missing columns: {sorted(missing)}"
            )

        df = df.merge(
            sacct[["jobid", "sacct_state", "sacct_exitcode"]],
            on="jobid",
            how="left",
        )

        ok = (df["sacct_state"] == "COMPLETED") & (df["sacct_exitcode"] == "0:0")
        if "rc" in df.columns:
            ok = ok & (pd.to_numeric(df["rc"], errors="coerce") == 0)

        df = df.loc[ok].copy()

    if return_counts:
        out = (
            df.groupby("method_with_flavor", dropna=False)
            .size()
            .reset_index(name="n_runs")
            .sort_values(["n_runs", "method_with_flavor"], ascending=[False, True])
            .reset_index(drop=True)
        )
        return out

    return sorted(df["method_with_flavor"].dropna().unique().tolist())
