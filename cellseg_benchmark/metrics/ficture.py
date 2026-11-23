#!/usr/bin/env python
"""FICTURE F1 metric: agreement between FICTURE factor identity and a segmentation's cell types."""

import gc
import logging
import pathlib

import joblib
import numpy as np
import pandas as pd
import scipy.spatial as ss
import seaborn as sns
import spatialdata as sd
from matplotlib import pyplot as plt
from spatialdata._io.io_shapes import _read_shapes

from .. import _constants
from .. import ficture_utils as fu
from .. import sdata_utils as su
from . import f1_score
from .utils import clean_method_name

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]: %(message)s")
logger = logging.getLogger("ficture_f1")


def compute_ficture_f1(
    adata,
    method,
    base_path,
    celltype_col="cell_type_revised",
    sample_col="sample",
    n_jobs=-1,
    **kwargs,
):
    """Per-cell-type F1 between FICTURE factors and a segmentation's cell types.

    Each transcript is labelled twice: with the cell type implied by its nearest
    FICTURE pixel's top factor (K1, within 5 um) and with the cell type of the
    segmentation cell it lies in. Both labels are mapped to the cell types in
    ``_constants.true_cluster`` before scoring, so labels such as Bergmann or
    Neurons-Dopa-Gaba count towards their parent, and transcripts with a label that
    is not scored (None or missing, e.g. Neurons-Other, Undefined) are dropped. Fits the
    standard ``compute_metric`` contract and scores all of a method's samples
    in parallel.

    Args:
        adata: Integrated AnnData for one method; needs ``sample_col`` and ``celltype_col`` in ``.obs``.
        method: Segmentation method name.
        base_path: Benchmark base directory.
        celltype_col: ``.obs`` column holding the segmentation cell types.
        sample_col: ``.obs`` column holding the sample names.
        n_jobs: Parallel workers over samples.
        **kwargs: Extra keyword arguments from the metric runner (ignored).

    Returns:
        DataFrame with one row per (sample, cell type) plus pooled ``all_samples`` rows.
        Columns: sample, cell_type, tp, fp, fn, precision, recall, f1. Besides one row
        per cell type, each sample block has four summary rows: f1_micro_all, f1_macro_all,
        f1_micro_vascular, f1_macro_vascular (micro = pooled counts, macro = mean of
        per-type f1; vascular = ECs/Pericytes/SMCs/VLMCs).
    """
    if any(method.startswith(p) for p in _constants.methods_3D):
        logger.info(f"[{method}] 3D method; skipping FICTURE F1.")
        return None

    obs = adata.obs[[sample_col, celltype_col]].copy()
    # obs ids are the boundary cell_id (first 10 chars) plus an AnnData suffix
    obs.index = obs.index.astype(str)
    if method in ["Baysor_2D_denovo", "Baysor_3D_denovo"]:
        obs.index = obs.index.str.extract(r"^(cell_\d+|[^_]+)", expand=False)
    else:
        if not obs.index.str.fullmatch(r"\d+").all():
            obs.index = obs.index.str.slice(0, 10)

    samples = obs[sample_col].unique().tolist()
    per_sample = joblib.Parallel(n_jobs=n_jobs, verbose=10)(
        joblib.delayed(_labelled_transcripts)(
            sample,
            obs.loc[obs[sample_col] == sample, celltype_col],
            method,
            base_path,
        )
        for sample in samples
    )

    labels = {sample: df for sample, df in per_sample if not df.empty}
    if not labels:
        logger.warning(f"No FICTURE F1 results for method '{method}'.")
        return pd.DataFrame()

    labels["all_samples"] = pd.concat(labels.values(), ignore_index=True)
    results = []
    for sample, df in labels.items():
        report = f1_score.compute_f1(df["ficture"], df["segmentation"]).reset_index()

        # micro (pooled counts) and macro (mean per-type f1), over all types and vascular only
        vascular = report[report["cell_type"].isin(_constants.vascular_celltypes)]
        summary = []
        for name, sub in [("all", report), ("vascular", vascular)]:
            tp, fp, fn = sub[["tp", "fp", "fn"]].sum()
            summary.append(
                {
                    "cell_type": f"f1_micro_{name}",
                    "tp": tp,
                    "fp": fp,
                    "fn": fn,
                    "precision": tp / (tp + fp) if tp + fp else 0.0,
                    "recall": tp / (tp + fn) if tp + fn else 0.0,
                    "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0,
                }
            )
            summary.append(
                {
                    "cell_type": f"f1_macro_{name}",
                    "f1": sub["f1"].mean() if len(sub) else 0.0,
                }
            )
        summary = pd.DataFrame(summary)

        report = pd.concat([report, summary], ignore_index=True)
        report.insert(0, "sample", sample)
        results.append(report)
    return pd.concat(results, ignore_index=True)


def _transcript_factors(sample, base_path):
    """Nearest FICTURE top-factor (K1, within 5 um) per transcript for one sample.

    Method-independent, so computed once and cached next to the FICTURE output as
    ``ficture_transcript_factors.parquet`` (columns: x, y, factor; factor -1 = no
    pixel within 5 um). Reused across all segmentation methods; recomputed if FICTURE
    was rerun after the cache was written.

    Raises:
        FileNotFoundError: If the sample has no FICTURE output.
    """
    cache = (
        pathlib.Path(base_path)
        / "samples"
        / sample
        / "results"
        / "Ficture"
        / "ficture_transcript_factors.parquet"
    )
    pixel_file = fu.find_ficture_output(sample, base_path)  # raises if missing
    if cache.exists() and cache.stat().st_mtime > pathlib.Path(pixel_file).stat().st_mtime:
        logger.info(f"[{sample}] loading cached transcript factors")
        return pd.read_parquet(cache)

    logger.info(f"[{sample}] computing transcript factors (no or outdated cache)")
    sdata = sd.read_zarr(
        pathlib.Path(base_path) / "samples" / sample / "sdata_z3.zarr",
        selection=("points",),
    )
    transcripts = sdata[f"{sample}_transcripts"][["x", "y"]].compute()
    del sdata

    meta = fu.parse_metadata(pixel_file)
    pixels = fu.read_ficture_pixels(pixel_file, usecols=["X", "Y", "K1"])
    xy = np.column_stack(
        [
            pixels["X"].to_numpy("float64") / float(meta["SCALE"]) + float(meta["OFFSET_X"]),
            pixels["Y"].to_numpy("float64") / float(meta["SCALE"]) + float(meta["OFFSET_Y"]),
        ]
    )
    k1 = pixels["K1"].to_numpy()
    del pixels
    distance, nearest = ss.cKDTree(xy).query(
        transcripts[["x", "y"]].to_numpy(), k=1, distance_upper_bound=5
    )
    del xy
    within = np.isfinite(distance)
    factor = np.full(len(transcripts), -1)
    factor[within] = k1[nearest[within]]

    # downcast to shrink the cache
    transcripts["x"] = transcripts["x"].astype("float32")
    transcripts["y"] = transcripts["y"].astype("float32")
    transcripts["factor"] = factor.astype("int8")
    transcripts.to_parquet(cache, compression="zstd", compression_level=3)
    return transcripts


def _labelled_transcripts(sample, celltypes, method, base_path):
    """Return per-transcript canonical labels (ficture, segmentation) for one sample."""
    try:
        transcripts = _transcript_factors(sample, base_path)
    except FileNotFoundError:
        logger.warning(f"[{sample}] no FICTURE output; skipping.")
        return sample, pd.DataFrame(columns=["ficture", "segmentation"])

    boundaries = sd.transform(
        _read_shapes(pathlib.Path(base_path) / "samples" / sample / "sdata_z3.zarr" / "shapes" / f"boundaries_{method}"),
        to_coordinate_system="micron",
    ).copy()

    ids = boundaries["cell_id"] if "cell_id" in boundaries.columns else boundaries.index
    boundaries["cell"] = pd.Index(ids).astype(str)

    # Observed label: segmentation cell type of the polygon each transcript lies in -> canonical.
    transcripts = su.assign_points_to_polygons(
        transcripts, boundaries, polygon_id_col="cell", output_col="cell"
    )
    del boundaries
    if not np.intersect1d(
        transcripts["cell"].dropna().unique(), celltypes.index.to_numpy()
    ).size:
        raise ValueError(
            f"[{sample}/{method}] no cell-id overlap between boundaries and obs"
        )
    labels = pd.DataFrame(
        {
            "ficture": transcripts["factor"]
            .map(fu.factor_cell_types(fu.find_ficture_output(sample, base_path)))
            .to_numpy(),
            "segmentation": transcripts["cell"]
            .map(dict(celltypes))
            .map(_constants.true_cluster)
            .to_numpy(),
        }
    ).dropna()

    del transcripts
    gc.collect()
    return sample, labels

def plot_ficture_f1(cohort, score="f1_macro_all", methods=None, show=True):
    """Dot plot of FICTURE F1 per method: mean (method color) over per-sample values (grey), best on top.

    Args:
        cohort: cohort name, or "merged" for all cohorts (mean weights cohorts equally).
        score: summary row to plot, one of f1_{macro,micro}_{all,vascular}.
        methods: optional list of clean method names to keep.
        show: show the figure.
    """
    cohorts = ["aging", "ABCAtlas", "VizgenMouseBrain"] if cohort == "merged" else [cohort]
    metrics = pathlib.Path(_constants.BASE_PATH) / "metrics"
    df = pd.concat([pd.read_csv(metrics / c / "ficture/ficture_f1.csv", index_col=0).assign(cohort=c) for c in cohorts])
    df = df[(df["sample"] != "all_samples") & (df.cell_type == score)]
    df["label"] = df.method.map(clean_method_name)
    if methods is not None:
        df = df[df.label.isin(methods)]
    means = df.groupby(["label", "cohort"]).f1.mean().groupby("label").mean().sort_values()
    y = pd.Series(range(len(means)), index=means.index)
    colors = df.drop_duplicates("label").set_index("label").method.map(lambda m: _constants.method_colors.get(m, "0.6"))
    df["y"] = df.label.map(y) + np.random.default_rng(0).uniform(-0.2, 0.2, len(df))

    fig, ax = plt.subplots(figsize=(4.5, 0.2 * len(means) + 1))
    for (c, d), marker in zip(df.groupby("cohort", sort=False), "os^"):
        ax.scatter(d.f1, d.y, s=7, marker=marker, color="0.7", lw=0, label=c)
    ax.scatter(means, y, s=40, color=colors[means.index], edgecolor="k", lw=0.5, zorder=3)
    ax.set_yticks(y, means.index, fontsize=8)
    ax.set_ylim(-0.7, len(means) - 0.3)
    ax.set(title=cohort, xlabel=f"FICTURE {score.split('_')[1]} F1" + (" (vascular)" if score.endswith("vascular") else ""))
    ax.grid(axis="x", lw=0.4, alpha=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    if cohort == "merged":
        ax.legend(frameon=False, fontsize=7, loc="lower right", handletextpad=0)
    sns.despine(ax=ax, left=True)
    (metrics / cohort / "ficture/plots").mkdir(parents=True, exist_ok=True)
    fig.savefig(metrics / cohort / f"ficture/plots/ficture_{score}.png", dpi=300, bbox_inches="tight")
    plt.show() if show else plt.close(fig)


def plot_ficture_f1_heatmap(cohort, show=True):
    """Heatmap of FICTURE F1 per method and cell type, samples pooled (mean over cohorts if merged)."""
    cohorts = ["aging", "ABCAtlas", "VizgenMouseBrain"] if cohort == "merged" else [cohort]
    metrics = pathlib.Path(_constants.BASE_PATH) / "metrics"
    df = pd.concat([pd.read_csv(metrics / c / "ficture/ficture_f1.csv", index_col=0) for c in cohorts])
    df = df[(df["sample"] == "all_samples") & ~df.cell_type.str.startswith("f1_")]
    mat = df.pivot_table(index=df.method.map(clean_method_name), columns="cell_type", values="f1")
    mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index, mat.mean().sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(0.4 * mat.shape[1] + 2.8, 0.24 * mat.shape[0] + 1.2))
    sns.heatmap(mat, annot=mat.apply(lambda col: col.map(lambda v: f"{v:.2f}".lstrip("0") if v == v else "")), fmt="",
                annot_kws={"size": 7.5}, cmap="YlOrRd", vmin=0, vmax=1, cbar_kws={"label": "F1", "shrink": 0.5}, ax=ax)
    ax.tick_params(labelsize=8, length=0)
    ax.set(xlabel="", ylabel="", title=f"{cohort}: FICTURE F1 per cell type (samples pooled)")
    (metrics / cohort / "ficture/plots").mkdir(parents=True, exist_ok=True)
    fig.savefig(metrics / cohort / "ficture/plots/ficture_f1_celltype_heatmap.png", dpi=300, bbox_inches="tight")
    plt.show() if show else plt.close(fig)
