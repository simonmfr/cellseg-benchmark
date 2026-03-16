import math
import os
import re
import subprocess
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
    base_path=BASE_PATH,
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
    if pass_method:
        results = metric_func(adata, method=method, base_path=base_path, **kwargs)
    else:
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


def chunked(lst, n=200):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

def normalize_jobname(jobname: str) -> str:
    if not isinstance(jobname, str):
        return ""
    # fix historic typo: vtp -> vpt
    return re.sub(r"^vtp", "vpt", jobname)

def parse_slurm_mem_to_gb(x: str) -> float:
    # sacct MaxRSS often: 1234K / 512M / 10G / 1T / 0 / ''
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
    mult = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4, "P": 1024**5}[unit]
    return (val * mult) / (1024**3)

def run_sacct(jobids):
    # Include ExitCode to filter success
    fmt = "JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,MaxRSS,JobName"
    cmd = ["sacct", "-X", "-P", "-n", "-j", ",".join(jobids), "--format", fmt]
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        raise RuntimeError(f"sacct failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    rows = []
    cols = fmt.split(",")
    for line in res.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("|")
        if len(parts) != len(cols):
            continue
        rows.append(dict(zip(cols, parts)))
    return pd.DataFrame(rows)

def method_with_flavor_from_row(jobname: str, key: str) -> str:
    """
    Create a single canonical string like:
      Baysor_CP1_DAPI_PolyT_0.8
      vpt2D_DAPI_PolyT
      Proseg_3D_vpt2_XYZ_vxl_3
      Proseg_3D_CP2_DAPI_PolyT_vxl_3
    """
    j = normalize_jobname(jobname)
    k = str(key)

    # --- Baysor_{key}_CP{cp}_{stain}_{conf} -> Baysor_CP{cp}_DAPI_{stain}_{conf}
    prefix = f"Baysor_{k}_"
    if j.startswith(prefix):
        rest = j[len(prefix):]  # CP1_PolyT_0.8
        m = re.match(r"^CP(?P<cp>\d+)_(?P<stain>[^_]+)_(?P<conf>.+)$", rest)
        if m:
            return f"Baysor_CP{m.group('cp')}_DAPI_{m.group('stain')}_{m.group('conf')}"
        return f"Baysor_{rest}"

    # --- vpt2D / vpt3D: vpt2D_{key}_{stain...} -> vpt2D_{stain...}
    for dim in ("2D", "3D"):
        prefix = f"vpt{dim}_{k}_"
        if j.startswith(prefix):
            rest = j[len(prefix):]  # DAPI_PolyT  OR  PolyT
            return f"vpt{dim}_{rest}"
        # sometimes only vpt2D_{key} (rare)
        if j == f"vpt{dim}_{k}":
            return f"vpt{dim}"

    # --- Proseg_3D vpt: Proseg_3D_{key}_vpt{dim}_{flavor}_vxl_{voxel}
    prefix = f"Proseg_3D_{k}_"
    if j.startswith(prefix) and "_vpt" in j:
        rest = j[len(prefix):]  # vpt2_flavor_vxl_...
        return f"Proseg_3D_{rest}"

    # --- Proseg_3D CP: Proseg_3D_{key}_CP{cp}_{stain}_vxl_{voxel} -> Proseg_3D_CP{cp}_DAPI_{stain}_vxl_{voxel}
    if j.startswith(prefix) and "_CP" in j:
        rest = j[len(prefix):]  # CP2_PolyT_vxl_3
        m = re.match(r"^CP(?P<cp>\d+)_(?P<stain>[^_]+)_vxl_(?P<voxel>.+)$", rest)
        if m:
            return f"Proseg_3D_CP{m.group('cp')}_DAPI_{m.group('stain')}_vxl_{m.group('voxel')}"
        return f"Proseg_3D_{rest}"

    # --- Proseg (2D) CP: Proseg_{key}_CP{cp}_{stain}_vxl_{voxel}
    prefix = f"Proseg_{k}_"
    if j.startswith(prefix) and "_CP" in j:
        rest = j[len(prefix):]  # CP2_PolyT_vxl_3
        m = re.match(r"^CP(?P<cp>\d+)_(?P<stain>[^_]+)_vxl_(?P<voxel>.+)$", rest)
        if m:
            return f"Proseg_CP{m.group('cp')}_DAPI_{m.group('stain')}_vxl_{m.group('voxel')}"
        return f"Proseg_{rest}"

    # --- Cellpose jobs: CP1_{key}_{stain} / CP2_{key}_{stain}
    prefix = f"CP1_{k}_"
    if j.startswith(prefix):
        stain = j[len(prefix):]
        return f"Cellpose_1_DAPI_{stain}"
    prefix = f"CP2_{k}_"
    if j.startswith(prefix):
        stain = j[len(prefix):]
        return f"Cellpose_2_DAPI_{stain}"

    # --- Simple methods (no extra flavor in jobname beyond key)
    simple = {
        f"voronoi_{k}": "voronoi",
        f"nuclei_{k}": "nuclei",
        f"merscope_{k}": "merscope_sdata",
        f"transcript_tif_{k}": "transcript_tif",
    }
    if j in simple:
        return simple[j]

    # --- rastered{width}_{key} -> rastered{width}
    m = re.match(rf"^rastered(?P<width>\d+)_({re.escape(k)})$", j)
    if m:
        return f"rastered{m.group('width')}"

    # fallback: keep normalized jobname (still useful)
    return j