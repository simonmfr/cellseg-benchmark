from pathlib import Path

import cellseg_benchmark as cb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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
    adata,
    obs_columns=None,
    obsm_columns=None,
    celltype_name="cell_type_revised",
    **kwargs,
):
    """Extract and save per-sample and per-celltype mean stats from adata.obs and obsm.

    Default behavior is to extract volume_final, area, sphericity, elongation, ovrlpy mean_integrity, PolyT and DAPI intensity

    Args:
        adata: anndata to extract morphology stats from
        obs_columns: names of obs columns to extract data from
        obsm_columns: dict with obsm names as keys and obsm column names as values.
        celltype_name: name of celltype column in adata.obs.
        kwargs: additional keyword arguments to pass.

    Returns:
        results DataFrame or None.
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
        Path(cb.BASE_PATH) / "metrics" / cohort / "general_stats" / "general_stats.csv"
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
        palette=cb._constants.method_colors,
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


def extract_mem_and_time(
    adata,
    method: str,
    ref_file_path: str | Path=Path(cb.BASE_PATH) / "misc/logs/job_runs.tsv",
    metrics_dir: str | Path=Path(cb.BASE_PATH) / "misc/extracted_job_stats",
    base_path=None,
    ignore_missing: bool=False,
    **kwargs,
) -> pd.DataFrame:
    """Read job metadata from ref_file_path and enrich it from the newest exported sacct TSV in metrics_dir.

    Args:
        adata: API compatibility.
        method (str): method name.
        ref_file_path (str or Path): path to reference TSV file with job information.
        metrics_dir (str or Path): path to metrics directory containing sacct read-outs.
        base_path: API compatibility.
        ignore_missing (bool): ignore methods without successful recorded segmentation.

    Returns:
        DataFrame with columns ["sample", "maxrss_gb", "elapsed_h", "alloccpus"]

    Notes:
        - Keeps only successful runs:
            sacct_state == COMPLETED
            sacct_exitcode == 0:0
            and rc == 0 if rc exists in the ref file
        - If the ref file contains repeated runs for the same sample+method,
          keeps the last successful one because the ref file is appended.
        - 'adata' and 'base_path' are unused and only kept for API compatibility.
    """

    def _missing_result(samples=None):
        """Build NaN dataframe for not successfully run or missing segmentations."""
        if samples is None:
            samples = pd.Series(dtype="object")
        else:
            samples = pd.Series(samples, dtype="object")
        return pd.DataFrame(
            {
                "sample": samples,
                "maxrss_gb": pd.Series([pd.NA] * len(samples), dtype="object"),
                "elapsed_h": pd.Series([pd.NA] * len(samples), dtype="object"),
                "alloccpus": pd.Series([pd.NA] * len(samples), dtype="object"),
            }
        ).reset_index(drop=True)

    ref = pd.read_csv(ref_file_path, sep="\t")
    ref["_ref_order"] = range(len(ref))

    ref["jobid"] = ref["jobid"].astype(str)
    ref["jobname"] = ref["jobname"].astype(str)
    ref["sample"] = ref["key"].astype(str)
    ref["jobname_norm"] = ref["jobname"].apply(cb.metrics.utils.normalize_jobname)

    ref["method_with_flavor"] = ref.apply(
        lambda r: cb.metrics.utils.method_with_flavor_from_row(r["jobname"], r["sample"]),
        axis=1,
    )

    ref = ref[ref["method_with_flavor"] == method].copy()
    if ref.empty:
        if ignore_missing:
            return _missing_result()
        raise LookupError(
            f"Method {method!r} not found in job file or not yet recorded."
        )

    latest_metrics_file = cb.metrics.utils.find_latest_job_data_tsv(metrics_dir)

    sacct = pd.read_csv(latest_metrics_file, sep="\t")
    if sacct.empty:
        raise ValueError(f"Latest metrics file is empty: {latest_metrics_file}")

    required_cols = {
        "jobid",
        "sacct_state",
        "sacct_exitcode",
        "elapsed_s",
        "alloccpus",
        "maxrss_gb",
    }
    missing = required_cols.difference(sacct.columns)
    if missing:
        raise ValueError(
            f"Metrics file {latest_metrics_file} is missing columns: {sorted(missing)}"
        )

    sacct["jobid"] = sacct["jobid"].astype(str)

    # avoid collisions with columns from ref file
    sacct = sacct[
        [
            "jobid",
            "sacct_state",
            "sacct_exitcode",
            "elapsed_s",
            "alloccpus",
            "maxrss_gb",
        ]
    ].rename(
        columns={
            "elapsed_s": "elapsed_s_sacct",
            "alloccpus": "alloccpus_sacct",
            "maxrss_gb": "maxrss_gb_sacct",
        }
    )

    df = ref.merge(sacct, on="jobid", how="left")

    ok = (df["sacct_state"] == "COMPLETED") & (df["sacct_exitcode"] == "0:0")
    if "rc" in df.columns:
        ok = ok & (pd.to_numeric(df["rc"], errors="coerce") == 0)

    df_ok = df.loc[ok].copy()
    if df_ok.empty:
        if ignore_missing:
            samples = (
                ref.sort_values("_ref_order")["sample"]
                .drop_duplicates()
                .tolist()
            )
            return _missing_result(samples)
        raise LookupError(
            f"No successful runs found for method {method!r} "
            f"in latest metrics file: {latest_metrics_file.name}"
        )

    # ref file is appended -> keep the last successful run per sample+method
    df_ok = (
        df_ok.sort_values("_ref_order")
        .drop_duplicates(subset=["sample", "method_with_flavor"], keep="last")
        .copy()
    )

    df_ok["elapsed_h"] = (
        pd.to_numeric(df_ok["elapsed_s_sacct"], errors="coerce") / 3600.0
    )

    out = df_ok[["sample", "maxrss_gb_sacct", "elapsed_h", "alloccpus_sacct"]].copy()

    out = out.rename(
        columns={
            "maxrss_gb_sacct": "maxrss_gb",
            "alloccpus_sacct": "alloccpus",
        }
    )

    out = out.reset_index(drop=True)
    return out

def plot_mem_and_time(cohort, metric, show: bool = False):
    """Violin plots of chosen metrics. Metrics can be "memory", "cpus", "duration"."""
    if metric not in ["memory", "cpus", "duration"]:
        raise ValueError(f"Metric {metric!r} is not supported. Chose one of memory, cpus or duration.")

    column_mapping = {
        "memory": "maxrss_gb",
        "cpus": "alloccpus",
        "duration": "elapsed_h",
    }
    col_name = column_mapping[metric]
    results_file = (
            Path(cb.BASE_PATH) / "metrics" / cohort / "Mem_and_time" / "mem_and_time.csv"
    )
    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    results_df = pd.read_csv(results_file, index_col=0)

    # Remove nan
    results_df = results_df[~results_df[col_name].isna()]

    # Remove outliers
    threshold = np.percentile(results_df[col_name], 99)
    results_df = results_df[results_df[col_name] <= threshold]

    dataset_order = results_df.groupby("method")[col_name].median().sort_values().index

    fig = plt.figure(figsize=(6, 6), dpi=300)
    plt.grid(True, alpha=0.3, zorder=0)
    sns.violinplot(
        results_df,
        y="method",
        x=col_name,
        hue="method",
        order=dataset_order,
        palette=cb._constants.method_colors,
        inner="quartile",
        linewidth=0.7,
        zorder=2,
        legend=False,
    )
    plt.tight_layout()
    if show:
        plt.show()
    fig.savefig(
        plot_path / f"{metric}.png", bbox_inches="tight"
    )
    plt.close(fig)
