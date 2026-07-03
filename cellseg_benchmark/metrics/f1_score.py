#!/usr/bin/env python
"""Per-label F1 from a confusion matrix over aligned true/predicted labels."""

import numpy as np
import pandas as pd


def compute_f1(true, pred, labels=None) -> pd.DataFrame:
    """Per-label precision, recall and F1 from aligned true/predicted labels.

    Rows where either label is missing (NaN or not in ``labels``) are ignored.

    Args:
        true: Iterable of ground-truth labels.
        pred: Iterable of predicted labels, aligned element-wise with ``true``.
        labels: Label order. Defaults to the sorted union of observed labels.

    Returns:
        DataFrame indexed by label (index name "cell_type") with columns
        tp, fp, fn, precision, recall, f1.
    """
    true = pd.Series(list(true), dtype="object")
    pred = pd.Series(list(pred), dtype="object")
    if labels is None:
        labels = sorted(set(true.dropna()) | set(pred.dropna()))
    labels = pd.Index(labels, name="cell_type")

    # int64: codes are int8 for <=127 categories, and codes * n would overflow.
    true_codes = pd.Categorical(true, categories=labels).codes.astype(np.int64)
    pred_codes = pd.Categorical(pred, categories=labels).codes.astype(np.int64)
    keep = (true_codes >= 0) & (pred_codes >= 0)
    n = len(labels)
    confusion = np.bincount(
        true_codes[keep] * n + pred_codes[keep], minlength=n * n
    ).reshape(n, n)

    tp = np.diag(confusion)
    fn = confusion.sum(axis=1) - tp
    fp = confusion.sum(axis=0) - tp
    with np.errstate(divide="ignore", invalid="ignore"):
        precision = np.where(tp + fp > 0, tp / (tp + fp), 0.0)
        recall = np.where(tp + fn > 0, tp / (tp + fn), 0.0)
        f1 = np.where(2 * tp + fp + fn > 0, 2 * tp / (2 * tp + fp + fn), 0.0)

    return pd.DataFrame(
        {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        },
        index=labels,
    )
