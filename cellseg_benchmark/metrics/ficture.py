#!/usr/bin/env python
"""FICTURE F1 metric: agreement between FICTURE factor identity and a segmentation's cell types."""

import gc
import logging
import pathlib

import joblib
import numpy as np
import pandas as pd
import scipy.spatial as ss
import spatialdata as sd

from .. import _constants
from .. import ficture_utils as fu
from .. import sdata_utils as su
from . import f1_score

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
    segmentation cell it lies in. Both labels are collapsed to canonical parent
    types (``_constants.true_cluster``) before scoring, so subtypes such as
    Astroependymal or Neurons-Dopa-Gaba count towards their parent, and transcripts
    touching markerless types (``_constants.unreliable_celltypes``) are dropped. Fits the
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

    factor_to_canonical = {
        int(factor): _constants.true_cluster[celltype]
        for factor, celltype in _constants.factor_to_celltype.items()
    }

    obs = adata.obs[[sample_col, celltype_col]].copy()
    # obs ids are the boundary cell_id (first 10 chars) plus an AnnData suffix
    obs.index = obs.index.astype(str)
    if not obs.index.str.fullmatch(r"\d+").all():
        obs.index = obs.index.str.slice(0, 10)

    samples = obs[sample_col].unique().tolist()
    per_sample = joblib.Parallel(n_jobs=n_jobs, verbose=10)(
        joblib.delayed(_labelled_transcripts)(
            sample,
            obs.loc[obs[sample_col] == sample, celltype_col],
            method,
            base_path,
            factor_to_canonical,
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
    pixel within 5 um). Reused across all segmentation methods.

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
    if cache.exists():
        logger.info(f"[{sample}] loading cached transcript factors")
        return pd.read_parquet(cache)

    logger.info(f"[{sample}] computing transcript factors (no cache)")
    pixel_file = fu.find_ficture_output(sample, base_path)  # raises if missing
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


def _labelled_transcripts(sample, celltypes, method, base_path, factor_to_canonical):
    """Return per-transcript canonical labels (ficture, segmentation) for one sample."""
    try:
        transcripts = _transcript_factors(sample, base_path)
    except FileNotFoundError:
        logger.warning(f"[{sample}] no FICTURE output; skipping.")
        return sample, pd.DataFrame(columns=["ficture", "segmentation"])

    sdata = sd.read_zarr(
        pathlib.Path(base_path) / "samples" / sample / "sdata_z3.zarr",
        selection=("shapes",),
    )
    boundaries = sd.transform(
        sdata[f"boundaries_{method}"], to_coordinate_system="micron"
    ).copy()
    del sdata

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
            "ficture": transcripts["factor"].map(factor_to_canonical).to_numpy(),
            "segmentation": transcripts["cell"]
            .map(dict(celltypes))
            .map(_constants.true_cluster)
            .to_numpy(),
        }
    )
    # drop transcripts touching markerless types (unreliable annotation) on either side
    labels = labels.replace(_constants.unreliable_celltypes, np.nan).dropna()

    del transcripts
    gc.collect()
    return sample, labels
