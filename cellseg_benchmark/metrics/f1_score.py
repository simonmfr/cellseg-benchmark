from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Literal, Union

import numpy as np
import pandas as pd
import spatialdata as sd
from joblib import Parallel, delayed
from numpy import ndarray
from scipy.spatial import cKDTree
from spatialdata import read_zarr

from cellseg_benchmark._constants import factor_to_celltype, true_cluster
from cellseg_benchmark.ficture_utils import (
    _find_ficture_output,
    _read_ficture_pixels,
    parse_metadata,
    process_coordinates,
)
from cellseg_benchmark.metrics.utils import _prepare_label_maps
from cellseg_benchmark.sdata_utils import _assign_points_to_polygons


def _f1_score(
    tp: Union[float, ndarray], fp: Union[float, ndarray], fn: Union[float, ndarray]
) -> Union[float, ndarray]:
    """Compute F1 score.

    Args:
        tp: True positive
        fp: False positive
        fn: False negative

    Returns:
        F1 score.
    """
    return (2 * tp) / (2 * tp + fp + fn)


def _compute_f1(
    data: pd.DataFrame | Dict[str, pd.DataFrame],
    flavor: Literal["f1", "micro", "macro", "all"] = "f1",
    celltype_col: str = "celltype",
    factor_col: str = "Factor",
    subset: List[str] | None = None,
    return_confusion: bool = False,
) -> pd.DataFrame:
    """Compute per-class F1 scores from a ground-truth celltype column and a
    predicted/assigned Factor column.

    Output is flipped:
    - columns = labels
    - rows = F1_score (and optionally TP, FP, FN, TN)

    Standard confusion matrix naming:
    - TP: true label == current class, predicted label == current class
    - FN: true label == current class, predicted label != current class
    - FP: true label != current class, predicted label == current class
    - TN: true label != current class, predicted label != current class.
    """
    if isinstance(data, pd.DataFrame):
        data = {"tmp": data}
    elif not isinstance(data, dict):
        raise TypeError("data must be a pandas DataFrame or a dict of DataFrames.")

    for key, df in data.items():
        if not isinstance(df, pd.DataFrame):
            raise TypeError(f"data['{key}'] is not a DataFrame.")
        missing = {celltype_col, factor_col} - set(df.columns)
        if missing:
            raise ValueError(
                f"data['{key}'] is missing required columns: {sorted(missing)}"
            )

    if subset is None:
        labels = set()
        for df in data.values():
            labels.update(df[celltype_col].dropna().unique())
            labels.update(df[factor_col].dropna().unique())
        labels = pd.Index(sorted(labels))
    else:
        labels = pd.Index(sorted(set(subset)))

    n = len(labels)
    cm_total = np.zeros((n, n), dtype=np.int64)

    for df in data.values():
        true = pd.Categorical(df[celltype_col], categories=labels)
        pred = pd.Categorical(df[factor_col], categories=labels)
        true_codes = true.codes.astype(np.int64)
        pred_codes = pred.codes.astype(np.int64)

        valid = (true_codes >= 0) & (pred_codes >= 0)
        cm = np.bincount(
            true_codes[valid] * n + pred_codes[valid],
            minlength=n * n,
        ).reshape(n, n)
        cm_total += cm

    tp = np.diag(cm_total)
    fn = cm_total.sum(axis=1) - tp
    fp = cm_total.sum(axis=0) - tp
    f1 = _f1_score(tp, fp, fn)
    out = pd.DataFrame([f1], index=["F1_score"], columns=labels)

    if return_confusion:
        total = cm_total.sum()
        tn = total - tp - fp - fn

        out.loc["TP"] = tp
        out.loc["FP"] = fp
        out.loc["FN"] = fn
        out.loc["TN"] = tn
        out = out.loc[["TP", "FP", "FN", "TN", "F1_score"]]

    if flavor in {"macro", "all"}:
        macro_f1 = float(np.mean(f1)) if len(f1) else 0.0
        out["macro F1_score"] = np.nan
        out.loc["F1_score", "macro F1_score"] = macro_f1

    if flavor in {"micro", "all"}:
        tp_total = int(tp.sum())
        fp_total = int(fp.sum())
        fn_total = int(fn.sum())
        denom = 2 * tp_total + fp_total + fn_total
        micro_f1 = 0.0 if denom == 0 else 2 * tp_total / denom

        out["micro F1_score"] = np.nan
        out.loc["F1_score", "micro F1_score"] = micro_f1

    return out


def f1_adapted(
    data: pd.DataFrame,
    min_value: Union[pd.DataFrame, float],
    correct_celltypes: Dict[str, str] = None,
) -> pd.DataFrame:
    """Computes the F1 score for mean and var information.

    Args:
        data: mean / var information.
        min_value: threshold, over which a cell is considered positive
        correct_celltypes: Optional. specify correct celltype per factor.

    Returns:
        F1 score per factor and over all factors to the given threshold.
    """
    if isinstance(min_value, float):
        min_value = pd.Series(min_value, index=data.columns)
    else:
        assert set(data.drop(columns=["cell type"]).columns) <= set(min_value.index), (
            "some cell types present in data are not represented in min_value DataFrame."
        )

    if correct_celltypes is None:
        f1_stats = pd.DataFrame(
            columns=data.drop(columns=["cell type"]).columns,
            index=data["cell type"].unique(),
        )

        for factor in data.drop(columns=["cell type"]).columns:
            for celltype in data["cell type"].unique():
                TP = (
                    data.loc[data["cell type"] == celltype, factor]
                    > float(min_value.loc[factor].values)
                ).sum()
                FP = (
                    data.loc[data["cell type"] != celltype, factor]
                    > float(min_value.loc[factor].values)
                ).sum()
                FN = len(data) - TP - FP
                f1_stats.loc[celltype, factor] = _f1_score(TP, FP, FN)
        f1_stats.sort_index(inplace=True)
        return f1_stats

    rates = {}

    for cluster, elements in correct_celltypes.items():
        rates[cluster] = defaultdict(float)
        rates[cluster]["TP"] = (
            data.loc[[x in elements for x in data["cell type"]], cluster]
            > min_value[cluster]
        ).sum()
        rates[cluster]["FP"] = (
            data.loc[[x not in elements for x in data["cell type"]], cluster]
            > min_value[cluster]
        ).sum()
        rates[cluster]["FN"] = len(data) - rates[cluster]["TP"] - rates[cluster]["FP"]
        rates[cluster]["F1_score"] = _f1_score(
            rates[cluster]["TP"], rates[cluster]["FP"], rates[cluster]["FN"]
        )

    f1_stats = defaultdict(float)
    for stats in rates.values():
        f1_stats["TP"] += stats["TP"]
        f1_stats["FP"] += stats["FP"]
        f1_stats["FN"] += stats["FN"]
    f1_stats["F1_score"] = _f1_score(f1_stats["TP"], f1_stats["FP"], f1_stats["FN"])

    stats_total = []
    names = []
    for cluster, stats in rates.items():
        names.append(cluster)
        stats_total.append(stats["F1_score"])
    names.append("Total")
    stats_total.append(f1_stats["F1_score"])
    stats_total = pd.DataFrame(stats_total, index=names, columns=["F1_statistics"]).T
    return stats_total


def _process_sample_ficture_f1(
    sample,
    obs_df,
    method,
    base_path,
    n_ficture,
    factor_to_celltype,
    true_cluster,
):
    """Compute ficture stats on one sample."""
    ficture_full_path = _find_ficture_output(sample, base_path, n_ficture)

    sdata = read_zarr(
        Path(base_path) / "samples" / sample / "sdata_z3.zarr",
        selection=("points", "shapes"),
    )

    points_data = sdata[f"{sample}_transcripts"].compute()
    boundaries = sd.transform(
        sdata[f"boundaries_{method}"],
        to_coordinate_system="micron",
    )

    del sdata
    boundaries = boundaries.copy()
    boundaries["p_id"] = boundaries.index

    ficture_pixels = _read_ficture_pixels(ficture_full_path)

    metadata = parse_metadata(ficture_full_path)
    ficture_pixels = process_coordinates(ficture_pixels, metadata)

    fic_microns = ficture_pixels[["x", "y"]].to_numpy()
    points_coords = points_data[["x", "y"]].to_numpy()

    tree = cKDTree(fic_microns)
    distances, indices = tree.query(points_coords, k=1, distance_upper_bound=5)

    result = points_data.copy()
    valid = np.isfinite(distances)
    result["assigned_factor"] = -1
    result["nearest_distance"] = distances
    result.loc[valid, "assigned_factor"] = ficture_pixels.iloc[indices[valid]][
        "K1"
    ].to_numpy()
    del ficture_pixels

    res = _assign_points_to_polygons(
        result,
        boundaries,
        polygon_id_col="p_id",
        output_col="assigned_polygon",
    )

    celltype = obs_df.loc[obs_df["sample"] == sample, "cell_type_revised"].copy()

    res_cts = res.merge(
        celltype,
        how="left",
        left_on="assigned_polygon",
        right_index=True,
    )

    factor_map, correct_celltypes = _prepare_label_maps(
        factor_to_celltype, true_cluster
    )

    res_cts["assigned_factor"] = res_cts["assigned_factor"].astype(str).map(factor_map)
    res_cts["cell_type_revised"] = res_cts["cell_type_revised"].map(correct_celltypes)
    eval_df = res_cts[["cell_type_revised", "assigned_factor"]].copy()

    metric = _compute_f1(
        eval_df,
        celltype_col="cell_type_revised",
        factor_col="assigned_factor",
        flavor="all",
        return_confusion=True,
    )

    sample_results = metric.T.reset_index().rename(
        columns={
            "index": "cell_type",
            "F1_score": "f1_score",
            "TP": "tp",
            "FP": "fp",
            "FN": "fn",
        }
    )
    sample_results = sample_results[
        ~sample_results["cell_type"].isin(["macro F1_score", "micro F1_score"])
    ].copy()

    sample_results["precision"] = sample_results["tp"] / (
        sample_results["tp"] + sample_results["fp"]
    )
    sample_results["recall"] = sample_results["tp"] / (
        sample_results["tp"] + sample_results["fn"]
    )

    sample_results["precision"] = sample_results["precision"].fillna(0)
    sample_results["recall"] = sample_results["recall"].fillna(0)
    sample_results.insert(0, "sample", sample)

    return sample, sample_results, eval_df


def ficture_f1_parallel(
    adata,
    method,
    base_path,
    n_ficture,
    n_jobs=-1,
    backend="loky",
):
    """Compute Ficture F1 score with parallelization."""
    # Do not mutate adata.obs_names in parallel code
    obs_df = adata.obs[["sample", "cell_type_revised"]].copy()
    obs_df.index = [x[:10] for x in adata.obs_names]

    samples = obs_df["sample"].unique().tolist()

    worker_out = Parallel(n_jobs=n_jobs, backend=backend, verbose=10)(
        delayed(_process_sample_ficture_f1)(
            sample=sample,
            obs_df=obs_df,
            method=method,
            base_path=base_path,
            n_ficture=n_ficture,
            factor_to_celltype=factor_to_celltype,
            true_cluster=true_cluster,
        )
        for sample in samples
    )

    per_sample_tables = []
    eval_frames = {}

    for sample, sample_results, eval_df in worker_out:
        per_sample_tables.append(sample_results)
        eval_frames[sample] = eval_df

    results = pd.concat(per_sample_tables, ignore_index=True)

    # Global score across all samples
    overall_metric = _compute_f1(
        eval_frames,
        celltype_col="cell_type_revised",
        factor_col="assigned_factor",
        flavor="all",
        return_confusion=True,
    )

    overall_results = overall_metric.T.reset_index().rename(
        columns={
            "index": "cell_type",
            "F1_score": "f1_score",
            "TP": "tp",
            "FP": "fp",
            "FN": "fn",
        }
    )
    overall_results = overall_results[
        ~overall_results["cell_type"].isin(["macro F1_score", "micro F1_score"])
    ].copy()

    overall_results["precision"] = overall_results["tp"] / (
        overall_results["tp"] + overall_results["fp"]
    )
    overall_results["recall"] = overall_results["tp"] / (
        overall_results["tp"] + overall_results["fn"]
    )

    overall_results["precision"] = overall_results["precision"].fillna(0)
    overall_results["recall"] = overall_results["recall"].fillna(0)
    overall_results.insert(0, "sample", "all_samples")

    results = pd.concat([results, overall_results], ignore_index=True)
    results.drop(results[results["cell_type"] == "unassigned"].index, inplace=True)

    return results
