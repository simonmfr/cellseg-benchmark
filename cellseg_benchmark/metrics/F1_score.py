import numpy as np
import pandas as pd
from collections import defaultdict


def compute_f1(data, correct_celltypes=None):
    """
    Computes the F1 score given a set of correct elements.
    Args:
        data: Data with factors and celltype annotation.
        correct_celltypes: specify correct celltype per factor.

    Returns:
        F1 score per factor and over all factors.
    """
    def f1_score(TP, FP, FN):
        return (2 * TP) / (2 * TP + FP + FN)

    if correct_celltypes is None:
        f1_stats = np.ones((len(data.drop(columns=["cell type"]).columns), len(data["cell type"].unique())))

        for i, factor in enumerate(data.drop(columns=["cell type"]).columns):
            for j, celltype in enumerate(data["cell type"].unique()):
                TP = data[factor][data["cell type"] == celltype].sum()
                FP = data[factor][data["cell type"] != celltype].sum()
                FN = 1 - TP - FP
                f1_stats[i, j] = f1_score(TP, FP, FN)

        f1_stats = pd.DataFrame(f1_stats, index=data["cell type"].unique(),
                                columns=data.drop(columns=["cell type"]).columns)
        f1_stats.sort_index(inplace=True)
        return f1_stats

    rates = {}

    for cluster, elements in correct_celltypes.items():
        rates[cluster] = defaultdict(float)
        rates[cluster]["TP"] = data[cluster][[x in elements for x in data["cell type"]]].sum()
        rates[cluster]["FP"] = data[cluster][[x not in elements for x in data["cell type"]]].sum()
        rates[cluster]["FN"] = 1-rates[cluster]["TP"]-rates[cluster]["FP"]
        rates[cluster]["F1_score"] = f1_score(rates[cluster]["TP"], rates[cluster]["FP"], rates[cluster]["FN"])

    f1_stats = defaultdict(float)
    for stats in rates.values():
        f1_stats["TP"] += stats["TP"]
        f1_stats["FP"] += stats["FP"]
        f1_stats["FN"] += stats["FN"]
    f1_stats["F1_score"] = f1_score(f1_stats["TP"], f1_stats["FP"], f1_stats["FN"])

    stats_total = []
    names = []
    for cluster, stats in rates.items():
        names.append(cluster)
        stats_total.append(stats["F1_score"])
    names.append("Total")
    stats_total.append(f1_stats["F1_score"])
    stats_total = pd.DataFrame(stats_total, index=names, columns=["F1_statistics"]).T
    return stats_total