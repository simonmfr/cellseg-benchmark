from typing import Dict, List, Literal, Union

import numpy as np
import pandas as pd
from numpy import ndarray


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


def compute_f1(
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
