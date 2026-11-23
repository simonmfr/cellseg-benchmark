from pathlib import Path

import anndata as ad
import geopandas as gpd
import h5py

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shapely

from . import utils
from .. import _constants
from ..adata_utils import plot_spatial_multiplot

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
    vascular_celltypes = ["ECs", "Pericytes", "SMCs", "VLMCs", "BAMs"]

    # compute for all cell types individually
    results = (
        df[["sample", celltype_name] + columns]
        .groupby(["sample", celltype_name], observed=True)
        .mean()
    )
    counts = (
        df[["sample", celltype_name] + columns]
        .groupby(["sample", celltype_name], observed=True)
        .size()
    )
    results = pd.concat([results, counts], axis=1)
    # compute for all cells together
    df_all = df[["sample"] + columns].groupby("sample", observed=True).mean()
    counts = df[["sample"] + columns].groupby("sample", observed=True).size()
    df_all = pd.concat([df_all, counts], axis=1)
    df_all[celltype_name] = "all"
    df_all = df_all.reset_index().set_index(["sample", celltype_name])
    results = pd.concat([results, df_all])
    # compute for only vascular subset
    df_vasc = df[df[celltype_name].isin(vascular_celltypes)][["sample"] + columns]
    counts = df_vasc.groupby("sample", observed=True).size()
    df_vasc = df_vasc.groupby("sample", observed=True).mean()
    df_vasc = pd.concat([df_vasc, counts], axis=1)
    df_vasc[celltype_name] = "vascular_subset"
    df_vasc = df_vasc.reset_index().set_index(["sample", celltype_name])
    results = pd.concat([results, df_vasc])
    results.rename({0: "n_cells"}, axis=1, inplace=True)
    return results.reset_index()


def tissue_polygons(cohort, base_path=_constants.BASE_PATH):
    """Cleaned BANKSY brain-region polygons merged per sample (manual outlier regions and stray islands removed)."""
    path = Path(base_path) / "misc" / "brain_regions" / f"{cohort}_brain_regions.parquet"
    return gpd.read_parquet(path).dissolve("sample").geometry


def _split_by_tissue(xy, tissue, buffer_um):
    """Yield (sample, cells, inside) with inside = cells within buffer_um of the tissue."""
    near = tissue.buffer(buffer_um)
    shapely.prepare(near.values)
    for s, d in xy[xy["sample"].isin(near.index)].groupby("sample"):
        yield s, d, shapely.contains_xy(near[s], d.x, d.y)


def compute_cell_density(adata, tissue, buffer_um=25, **kwargs):
    """Post-QC cells inside the tissue, tissue area (mm²) and cells per mm², per sample and for all samples.

    Cells more than buffer_um outside the tissue polygons, e.g. debris around the section, are not counted.
    """
    xy = pd.DataFrame(adata.obsm["spatial"][:, :2], columns=["x", "y"]).assign(
        sample=adata.obs["sample"].astype(str).to_numpy()
    )
    df = pd.DataFrame([{"sample": s, "n_cells": inside.sum(), "tissue_mm2": tissue[s].area / 1e6}
                       for s, _, inside in _split_by_tissue(xy, tissue, buffer_um)])
    df = pd.concat([df, df[["n_cells", "tissue_mm2"]].sum().to_frame().T.assign(sample="all")], ignore_index=True)
    return df.astype({"n_cells": int}).assign(cells_per_mm2=df.n_cells / df.tissue_mm2)


def plot_cell_density(cohort, tissue, path, buffer_um=25, base_path=_constants.BASE_PATH):
    """QC plot per sample: tissue border, kept cells (subsampled) and dropped cells of all methods."""
    pts = [pd.DataFrame(shapely.get_coordinates(shapely.segmentize(g.boundary, 10)), columns=["x", "y"])
           .assign(sample=s, status="border") for s, g in tissue.items()]
    for f in (Path(base_path) / "analysis" / cohort).glob("*/adatas/adata_integrated.h5ad.gz"):
        with h5py.File(f, "r") as h:
            s = h["obs/sample"]
            sample = s["categories"].asstr()[:][s["codes"][:]] if isinstance(s, h5py.Group) else s.asstr()[:]
            xy = pd.DataFrame(h["obsm/spatial"][:, :2], columns=["x", "y"]).assign(sample=sample)
        for _, d, inside in _split_by_tissue(xy, tissue, buffer_um):
            pts += [d[~inside].assign(status="dropped"),
                    d[inside].sample(min(2000, inside.sum()), random_state=0).assign(status="in tissue")]
    pts = pd.concat(pts, ignore_index=True)
    qc = ad.AnnData(obs=pts[["sample", "status"]].astype("category").set_index(pts.index.astype(str)))
    qc.obsm["spatial"] = pts[["x", "y"]].to_numpy()
    path = Path(path)
    plot_spatial_multiplot(qc, "status", path.parent, save_name=path.name, sort=True,
                           palette={"in tissue": "lightgrey", "border": "black", "dropped": "red"},
                           max_points_per_sample=len(pts))


def extract_general_stats(
    adata,
    obs_columns=None,
    obsm_columns=None,
    celltype_name="cell_type_revised",
    **kwargs,
):
    """Extract and save per-sample and per-celltype mean stats from adata.obs and obsm.

    Default behavior is to extract volume_final, area, circularity, elongation, n_counts, PolyT and DAPI intensity

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
        obs_columns = ["volume_final", "area", "circularity", "sphericity_3d", "elongation", "n_counts"]
    if obsm_columns is None:
        obsm_columns = {"intensities": ["PolyT", "DAPI"]}
    # prepare adata by putting obsm columns in obs
    df = adata.obs.copy()
    df = df.reindex(columns=df.columns.union(obs_columns, sort=False))
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
        Path(_constants.BASE_PATH) / "metrics" / cohort / "general_stats" / "general_stats.csv"
    )
    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    results_df = pd.read_csv(results_file, index_col=0)
    # select those with selected celltype
    results_df = results_df[results_df["cell_type_revised"] == celltype]
    results_df['method'] = results_df['method'].map(utils.clean_method_name)

    palette = {utils.clean_method_name(key): value for key, value in _constants.method_colors.items()}

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
        palette=palette,
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
    ref_file_path: str | Path=Path(_constants.BASE_PATH) / "misc/logs/run_log.tsv",
    metrics_dir: str | Path=Path(_constants.BASE_PATH) / "misc/extracted_job_stats",
    base_path=None,
    ignore_missing: bool=False,
    **kwargs,
) -> pd.DataFrame:
    """Runtime and peak memory per sample of one method, from run_log.tsv and the exported sacct tables.

    Args:
        adata: used to find the cohort.
        method (str): method name, matched to the results folder of each logged run.
        ref_file_path (str or Path): run_log.tsv.
        metrics_dir (str or Path): folder with the *_job_data.tsv exports of export_job_metrics_tsv.
        base_path: API compatibility.
        ignore_missing (bool): return an empty table for methods without logged runs.

    Returns:
        DataFrame with columns ["sample", "maxrss_gb", "elapsed_h", "alloccpus"]

    Notes:
        - Test runs (*_keep_cache) are left out. CoolMUC-4 job IDs clash with the AI system; their sacct values come from cm4_job_data.tsv, which sorts after the dated exports and so wins. Proseg_3D_vpt3D_* folders were later renamed Proseg_3D_vpt3D_DAPI_*. Per sample, the last COMPLETED run counts.
        - VPT jobs failed on update-vzg after writing all outputs until 48dffb6, so a FAILED last VPT run counts as done.
        - ComSeg resumes from cached patches, and since 2026-10-04 timed-out VPT runs were finished by
          recovery jobs on the missing tiles, so runs that did not complete directly before the counted
          run are added to its runtime (memory and CPUs: max).
    """
    cohorts = {s.split("_")[0] for s in adata.obs["sample"].unique()}
    assert len(cohorts) == 1, "more than one cohort found. Cohort recognition is sensitive to '_'"
    ref = pd.read_csv(ref_file_path, sep="\t")
    ref = ref[ref["key"].str.startswith(cohorts.pop())
              & (ref["result_dir"].str.extract(r"/results/([^/]+)", expand=False)
                 .str.replace(r"^Proseg_3D_vpt3D_(?!DAPI_)", "Proseg_3D_vpt3D_DAPI_", regex=True) == method)
              & ~ref["jobname"].str.endswith("_keep_cache")]
    if ref.empty:
        if ignore_missing:
            return pd.DataFrame(columns=["sample", "maxrss_gb", "elapsed_h", "alloccpus"])
        raise LookupError(f"Method {method!r} not found in job file or not yet recorded.")

    sacct = pd.concat(pd.read_csv(p, sep="\t") for p in sorted(Path(metrics_dir).glob("*_job_data.tsv")))
    runs = ref.merge(sacct.groupby("jobid", as_index=False).last(), on="jobid", suffixes=("_log", "")).sort_values("start_iso")
    since = "" if method.startswith("ComSeg") else "2026-10-04" if method.startswith("vpt_") else "9999"

    def counted(g):
        ok = (g["sacct_state"] == "COMPLETED").to_numpy()
        ok[-1] |= method.startswith("vpt_") and g["sacct_state"].iloc[-1] == "FAILED"
        done = np.flatnonzero(ok)
        if not len(done):
            return None
        g = g.iloc[(done[-2] + 1 if len(done) > 1 else 0):done[-1] + 1]
        g = g[(g["start_iso"] >= since) | (np.arange(len(g)) == len(g) - 1)]
        return {"sample": g["key"].iloc[0], "maxrss_gb": g["maxrss_gb"].max(), "elapsed_h": g["elapsed_s"].sum() / 3600, "alloccpus": g["alloccpus"].max()}

    rows = [r for _, g in runs.groupby("key") if (r := counted(g))]
    return pd.DataFrame(rows, columns=["sample", "maxrss_gb", "elapsed_h", "alloccpus"])


def plot_mem_and_time(cohort, metric=None, show: bool = False):
    """Violin plots of chosen metrics. Metrics can be "memory", "cpus", "duration"."""
    if isinstance(metric, str):
        if metric not in ["memory", "cpus", "duration"]:
            raise ValueError(f"Metric {metric!r} is not supported. Choose one of memory, cpus or duration.")
        metric = [metric]
    elif isinstance(metric, list):
        if not all([x in ["memory", "cpus", "duration"] for x in metric]):
            raise ValueError(f"Metric {metric!r} is not supported. Choose subset of memory, cpus or duration.")
    if metric is None:
        metric = ["memory", "cpus", "duration"]

    column_mapping = {
        "memory": "maxrss_gb",
        "cpus": "alloccpus",
        "duration": "elapsed_h",
    }

    for m in metric:
        col_name = column_mapping[m]
        results_file = (
                Path(_constants.BASE_PATH) / "metrics" / cohort / "Mem_and_time" / "mem_and_time.csv"
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
        palette = {utils.clean_method_name(m): _constants.method_colors[m] for m in dataset_order}
        dataset_order = [utils.clean_method_name(m) for m in dataset_order]

        results_df['method'] = results_df['method'].map(utils.clean_method_name)

        fig = plt.figure(figsize=(6, 6), dpi=300)
        plt.grid(True, alpha=0.3, zorder=0)
        sns.violinplot(
            results_df,
            y="method",
            x=col_name,
            hue="method",
            order=dataset_order,
            palette=palette,
            inner="quartile",
            linewidth=0.7,
            zorder=2,
            legend=False,
        )
        plt.tight_layout()
        if show:
            plt.show()
        fig.savefig(
            plot_path / f"{m}.png", bbox_inches="tight"
        )
        plt.close(fig)
