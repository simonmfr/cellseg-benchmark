from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from cellseg_benchmark import BASE_PATH
from cellseg_benchmark._constants import method_colors
from cellseg_benchmark.metrics.utils import (
    chunked,
    normalize_jobname,
    parse_slurm_mem_to_gb,
    run_sacct,
    method_with_flavor_from_row
)


def _extract_stats(df, columns, celltype_name="cell_type_revised"):
    """Extract and save per-sample and per-celltype mean stats from adata.obs.

    Args:
        df: dataframe to extract morphology stats from
        celltype_name: name of celltype column in df.
        columns: names of columns to extract data from
    Returns:
        results DataFrame or None
    """
    if celltype_name not in df.columns:
        print(f"{celltype_name} not found in adata.obs")
        return None

    # celltypes to use for vascular subset
    vascular_celltypes = ["ECs", "Pericytes", "SMCs", "VLMCs"]

    # compute for all cell types individually
    results = (
        df[["sample", celltype_name] + columns]
        .groupby(["sample", celltype_name], observed=True)
        .mean()
    )
    # compute for all cells together
    df_all = df[["sample"] + columns].groupby("sample", observed=True).mean()
    df_all[celltype_name] = "all"
    df_all = df_all.reset_index().set_index(["sample", celltype_name])
    results = pd.concat([results, df_all])
    # compute for only vascular subset
    df_vasc = df[df[celltype_name].isin(vascular_celltypes)][["sample"] + columns]
    df_vasc = df_vasc.groupby("sample", observed=True).mean()
    df_vasc[celltype_name] = "vascular_subset"
    df_vasc = df_vasc.reset_index().set_index(["sample", celltype_name])
    results = pd.concat([results, df_vasc])
    return results.reset_index()


def extract_general_stats(
    adata, obs_columns=None, obsm_columns=None, celltype_name="cell_type_revised", **kwargs
):
    """Extract and save per-sample and per-celltype mean stats from adata.obs and obsm.

    Default behavior is to extract volume_final, area, sphericity, elongation, ovrlpy mean_integrity, PolyT and DAPI intensity

    Args:
        adata: anndata to extract morphology stats from
        obs_columns: names of obs columns to extract data from
        obsm_columns: dict with obsm names as keys and obsm column names as values.
        celltype_name: name of celltype column in adata.obs.

    Returns:
        results DataFrame or None
    """
    # set default values
    if obs_columns is None:
        obs_columns = ["volume_final", "area", "sphericity", "elongation"]
    if obsm_columns is None:
        obsm_columns = {
            "intensities": ["PolyT", "DAPI"],
            "Ovrlpy_stats": ["mean_integrity"],
        }
    # prepare adata by putting obsm columns in obs
    df = adata.obs.copy()
    for key, values in obsm_columns.items():
        for value in values:
            new_key = f"{key}_{value}"
            try:
                df[new_key] = adata.obsm[key][value]
                obs_columns.append(new_key)
            except KeyError:
                print(f"{key}/{value} not found in obsm! Skipping {key}/{value}.")
                continue
    return _extract_stats(df, obs_columns, celltype_name)


def plot_general_stats(cohort, metric, celltype="all", show=False):
    """Plot general stats."""
    results_file = (
        Path(BASE_PATH) / "metrics" / cohort / "general_stats" / "general_stats.csv"
    )
    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    results_df = pd.read_csv(results_file, index_col=0)
    # select those with selected celltype
    results_df = results_df[results_df["cell_type_revised"] == celltype]

    # Remove nan
    results_df = results_df[~results_df[metric].isna()]

    # Remove outliers
    threshold = np.percentile(results_df[metric], 99)
    results_df = results_df[results_df[metric] <= threshold]

    dataset_order = results_df.groupby("method")[metric].median().sort_values().index

    fig = plt.figure(figsize=(6, 6), dpi=300)
    plt.grid(True, alpha=0.3, zorder=0)
    sns.violinplot(
        results_df,
        y="method",
        x=metric,
        hue="method",
        order=dataset_order,
        palette=method_colors,
        inner="quartile",
        linewidth=0.7,
        zorder=2,
        legend=False,
    )
    plt.tight_layout()
    if show:
        plt.show()
    fig.savefig(
        plot_path / f"general_stats_{metric}_{celltype}.png", bbox_inches="tight"
    )
    plt.close(fig)


def extract_mem_and_time(adata, ref_file_path, method):
    df = pd.read_csv(ref_file_path, sep="\t")
    df["jobid"] = df["jobid"].astype(str)
    df["jobname"] = df["jobname"].astype(str)
    df["jobname_norm"] = df["jobname"].apply(normalize_jobname)
    df["sample"] = df["key"].astype(str)

    df["method_with_flavor"] = df.apply(
        lambda r: method_with_flavor_from_row(r["jobname_norm"], r["sample"]),
        axis=1
    )
    df = df[df["method_with_flavor"] == method].copy()
    if df.empty:
        raise LookupError(f"Method {method} not found in job file or not yet recorded.")

    # method_with_flavor column
    df["method_with_flavor"] = df.apply(lambda r: method_with_flavor_from_row(r["jobname_norm"], r["sample"]), axis=1)

    # ---------- sacct enrichment ----------
    jobids = sorted(df["jobid"].dropna().unique().tolist())
    sacct_frames = []
    for ch in chunked(jobids, 200):
        sacct_frames.append(run_sacct(ch))
    sacct = pd.concat(sacct_frames, ignore_index=True)

    # normalize jobid base (strip steps: 12345.batch -> 12345)
    sacct["jobid_base"] = sacct["JobIDRaw"].astype(str).str.split(".").str[0]

    # parse metrics
    sacct["Elapsed_s"] = pd.to_numeric(sacct["ElapsedRaw"], errors="coerce")
    sacct["AllocCPUS_n"] = pd.to_numeric(sacct["AllocCPUS"], errors="coerce")
    sacct["MaxRSS_GB"] = sacct["MaxRSS"].apply(parse_slurm_mem_to_gb)

    # prefer MaxRSS from .batch if available, but keep elapsed/cpus from max over steps
    is_batch = sacct["JobIDRaw"].astype(str).str.endswith(".batch")

    agg_all = sacct.groupby("jobid_base", as_index=False).agg(
        sacct_state=("State", "first"),
        sacct_exitcode=("ExitCode", "first"),
        elapsed_s=("Elapsed_s", "max"),
        alloccpus=("AllocCPUS_n", "max"),
    )

    agg_batch = sacct[is_batch].groupby("jobid_base", as_index=False).agg(
        maxrss_gb=("MaxRSS_GB", "max"),
    )

    agg = agg_all.merge(agg_batch, on="jobid_base", how="left")

    df = df.merge(agg, left_on="jobid", right_on="jobid_base", how="left").drop(columns=["jobid_base"])

    # ---------- keep only successful segmentations ----------
    # success: sacct says COMPLETED + ExitCode 0:0 (and optionally your script rc == 0 if present)
    ok = (df["sacct_state"] == "COMPLETED") & (df["sacct_exitcode"] == "0:0")
    if "rc" in df.columns:
        ok = ok & (pd.to_numeric(df["rc"], errors="coerce") == 0)

    df_ok = df.loc[ok].copy()
    df_ok["elapsed_h"] = df_ok["elapsed_s"] / 3600.0
    out = df_ok[["sample", "method_with_flavor", "maxrss_gb", "elapsed_h", "alloccpus"]].copy()
    out = out.reset_index().set_index(["sample"])
    return out
