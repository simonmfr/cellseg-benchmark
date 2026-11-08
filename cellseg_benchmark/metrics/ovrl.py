"""Z integrity: share of a segmentation's signal that ovrlpy flags as mixed across z."""

import logging
import os
import pathlib

import numpy as np
import ovrlpy
import pandas as pd
import polars
import spatialdata as sd
import yaml
from rasterio.features import rasterize

from .. import _constants

NPZ = "vertical_doublets_ovrlpy_output.npz"
OFFSET = 13


def compute_ovrl(
    sample: str, sample_dir: str, data_dir: str, logger: logging.Logger = None, **kwargs
) -> None:
    """Compute ovrlpy output.

    Args:
        sample: name of sample
        sample_dir: Sample directory
        data_dir: Directory to transcript .csv.
        logger: logging.Logger instance.

    Returns:
        None
    """
    if not os.path.exists(
        os.path.join(sample_dir, NPZ)
    ):
        coords_df = pd.read_csv(
            os.path.join(data_dir, "detected_transcripts.csv"), index_col=0
        )
        if "global_x" in coords_df.columns and "global_y" in coords_df.columns:
            if logger is not None:
                logger.info("Basing ovrlpy on global_x, global_y and global_z columns.")
            coords_df = coords_df[["gene", "global_x", "global_y", "global_z"]].rename(
                columns={"global_x": "x", "global_y": "y", "global_z": "z"}
            )
        else:
            if logger is not None:
                logger.info("Basing ovrlpy on x, y and global_z columns.")
            coords_df = coords_df[["gene", "x", "y", "global_z"]].rename(columns={"global_z": "z"})
        run_ovrlpy(sample, coords_df, sample_dir)
    else:
        if logger is not None:
            logger.warning(f"Ovrlpy output exists, skipping ovrlpy run for {sample}")
        else:
            print("Ovrlpy output exists, skipping ovrlpy run.")


def run_ovrlpy(
    sample_name: str,
    coordinate_df: pd.DataFrame | polars.DataFrame,
    data_dir: str,
    logger: logging.Logger = None,
) -> None:
    """Processes MERSCOPE sample using ovrlpy to detect vertical doublets.

    Reads transcript coordinates, runs ovrlpy analysis,
    and saves the integrity and signal matrices as a compressed .npz file.

    Args:
        sample_name: Name of the MERSCOPE sample.
        coordinate_df: Coordinate DataFrame of transcript coordinates.
        data_dir: Base directory to save output files.
        logger: logging.Logger instance.

    Returns:
        None
    """
    n_workers = os.cpu_count()

    if logger:
        logger.info(f"[{sample_name}] Starting analysis...")
    else:
        print(f"[{sample_name}] Starting analysis...")

    ovrlpy_obj = ovrlpy.Ovrlp(coordinate_df, n_workers=n_workers, random_state=42)
    ovrlpy_obj.analyse()

    output_path = os.path.join(data_dir, NPZ)
    np.savez_compressed(
        output_path, integrity=ovrlpy_obj.integrity_map, signal=ovrlpy_obj.signal_map
    )

    if logger:
        logger.info(
            f"[{sample_name}] Analysis complete. Results saved to: data_dir/vertical_doublets_ovrlpy_output.npz"
        )
    else:
        print(
            f"[{sample_name}] Analysis complete. Results saved to: data_dir/vertical_doublets_ovrlpy_output.npz"
        )


def _vsi_maps(sample, base_path):
    """VSI and signal maps of one sample, plus the micron to VSI pixel shift."""
    npz = np.load(pathlib.Path(base_path) / "samples" / sample / NPZ)
    meta = yaml.safe_load(open(pathlib.Path(base_path) / "misc" / "sample_metadata.yaml"))
    tf = np.loadtxt(pathlib.Path(meta[sample]["path"], "images", "micron_to_mosaic_pixel_transform.csv"))
    shift = tf[[0, 1], 2] / tf[[0, 1], [0, 1]] - OFFSET
    return npz["integrity"], npz["signal"], shift


def _cells(sample, method, cell_ids, shift, base_path):
    """Outlines of the given cells in VSI pixel coordinates, indexed by cell id."""
    sdata = sd.read_zarr(pathlib.Path(base_path) / "samples" / sample / "sdata_z3.zarr", selection=("shapes",))
    cells = sd.transform(sdata[f"boundaries_{method}"], to_coordinate_system="micron")
    if "radius" in cells:
        cells.geometry = cells.geometry.buffer(cells["radius"])
    geoms = cells.geometry.translate(*shift)
    geoms.index = pd.Index((cells["cell_id"] if "cell_id" in cells else cells.index.to_series()).astype(str))
    cell_ids = pd.Index(cell_ids).astype(str)
    stripped = cell_ids.str.replace(r"-\d+$", "", regex=True)
    key = pd.Index(np.where(cell_ids.isin(geoms.index), cell_ids, stripped))
    to_obs = dict(zip(key, cell_ids))
    out = geoms[geoms.index.isin(key) & ~(geoms.is_empty | geoms.isna())]
    out.index = out.index.map(to_obs)
    return out


def compute_z_integrity(
    adata,
    method,
    base_path,
    sample_col="sample",
    signal_thresholds=(1, 2, 3),
    vsi_thresholds=(0.6, 0.7, 0.8),
    **kwargs,
):
    """Share of low-VSI signal inside a method's cells, per sample.

    Pools all pixels inside the outlines of the cells in ``adata`` (the union of
    footprints), keeps pixels with ovrlpy signal above the signal threshold, and
    reports the percentage with VSI below the VSI threshold. Pooling makes the
    score independent of how a region is split into cells; lower is better.
    Coverage is the share of all signal pixels inside the outlines. Fits the
    ``compute_metric`` contract; 3D methods are skipped.

    Args:
        adata: Integrated AnnData for one method; needs ``sample_col`` in ``.obs``.
        method: Segmentation method name.
        base_path: Benchmark base directory.
        sample_col: ``.obs`` column holding the sample names.
        signal_thresholds: Minimum ovrlpy signal of a pixel to count.
        vsi_thresholds: VSI below which a pixel counts as low.
        **kwargs: Extra keyword arguments from the metric runner (ignored).

    Returns:
        DataFrame with one row per sample and threshold pair, plus ``all_samples``
        rows with the mean over samples. Columns: sample, signal_threshold,
        vsi_threshold, n_cells, coverage_pct, mean_vsi, low_vsi_pct.
    """
    if any(method.startswith(p) for p in _constants.methods_3D):
        return None
    rows = []
    for sample, obs in adata.obs.groupby(sample_col, observed=True):
        vsi, signal, shift = _vsi_maps(sample, base_path)
        cells = _cells(sample, method, obs.index, shift, base_path)
        inside = rasterize(((g, 1) for g in cells), out_shape=vsi.shape, dtype=np.uint8) > 0
        for s in signal_thresholds:
            tissue = signal > s
            vals = vsi[inside & tissue]
            for v in vsi_thresholds:
                rows.append({"sample": sample, "signal_threshold": s, "vsi_threshold": v, "n_cells": cells.index.nunique(),
                             "coverage_pct": 100 * vals.size / tissue.sum(), "mean_vsi": vals.mean() if vals.size else np.nan,
                             "low_vsi_pct": 100 * (vals < v).mean() if vals.size else np.nan})
    df = pd.DataFrame(rows)
    mean = df.groupby(["signal_threshold", "vsi_threshold"], as_index=False).mean(numeric_only=True).assign(sample="all_samples")
    return pd.concat([df, mean[df.columns]], ignore_index=True)


def low_vsi_per_cell(adata, method, base_path, sample_col="sample", signal_threshold=2, vsi_threshold=0.7):
    """Per cell, the share of its signal pixels with low VSI.

    Used to check that low-VSI cells carry mixed signal (e.g. lower marker F1).
    Where outlines overlap, a pixel counts for the cell drawn last.

    Returns:
        Series indexed like ``adata.obs`` (NaN for cells without signal pixels).
    """
    out = []
    for sample, obs in adata.obs.groupby(sample_col, observed=True):
        vsi, signal, shift = _vsi_maps(sample, base_path)
        cells = _cells(sample, method, obs.index, shift, base_path)
        labels = rasterize(zip(cells, range(1, len(cells) + 1)), out_shape=vsi.shape, dtype=np.int32)
        on = (labels > 0) & (signal > signal_threshold)
        n = np.bincount(labels[on], minlength=len(cells) + 1)[1:]
        low = np.bincount(labels[on], weights=vsi[on] < vsi_threshold, minlength=len(cells) + 1)[1:]
        counts = pd.DataFrame({"low": low, "n": n}, index=cells.index).groupby(level=0).sum()
        out.append(counts["low"] / counts["n"].replace(0, np.nan))
    return pd.concat(out).reindex(adata.obs.index.astype(str))


def low_vsi_per_gene(sample, base_path, signal_threshold=2, vsi_threshold=0.7):
    """Per gene, the number of transcripts on signal pixels and how many fall on low-VSI pixels.

    Segmentation free: uses only the transcripts and the ovrlpy maps of one sample.

    Returns:
        DataFrame with columns gene, n, n_low.
    """
    vsi, signal, shift = _vsi_maps(sample, base_path)
    meta = yaml.safe_load(open(pathlib.Path(base_path) / "misc" / "sample_metadata.yaml"))
    csv = pathlib.Path(meta[sample]["path"], "detected_transcripts.csv")
    xy = ["global_x", "global_y"] if "global_x" in pd.read_csv(csv, nrows=0).columns else ["x", "y"]
    parts = []
    for tx in pd.read_csv(csv, usecols=[*xy, "gene"], dtype={xy[0]: "float32", xy[1]: "float32", "gene": "category"}, chunksize=5_000_000):
        c = np.floor(tx[xy[0]].to_numpy() + shift[0]).astype(int)
        r = np.floor(tx[xy[1]].to_numpy() + shift[1]).astype(int)
        ok = (r >= 0) & (r < vsi.shape[0]) & (c >= 0) & (c < vsi.shape[1])
        r, c, gene = r[ok], c[ok], tx["gene"].to_numpy()[ok]
        on = signal[r, c] > signal_threshold
        df = pd.DataFrame({"gene": np.asarray(gene)[on], "low": vsi[r[on], c[on]] < vsi_threshold})
        parts.append(df.groupby("gene", observed=True)["low"].agg(n="size", n_low="sum"))
    return pd.concat(parts).groupby(level=0).sum().rename_axis("gene").reset_index()
