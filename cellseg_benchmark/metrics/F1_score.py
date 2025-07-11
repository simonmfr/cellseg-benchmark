from typing import Union, Literal, List, Dict

from numpy import ndarray
import pandas as pd
from collections import defaultdict
from statistics import mean

def _f1_score(tp: Union[float, ndarray], fp: Union[float, ndarray], fn: Union[float, ndarray]) -> Union[float, ndarray]:
    return (2 * tp) / (2 * tp + fp + fn)

def compute_f1(data: pd.DataFrame, general_stats: pd.DataFrame=None,
               flavor: Literal["f1", "micro", "macro", "all"]="f1",
               correct_celltypes: Dict[str, List[str]]=None, weighted: bool=False, subset: List[str] = None
               ) -> pd.DataFrame:
    """
    Computes the F1 score given area information.
    Args:
        data: Data with factors and celltype annotation.
        general_stats: If provided, data must be total for FN computation. Otherwise, data is relative
        correct_celltypes: specify correct celltype per factor.
        weighted: Data based on weighted or unweighted area.

    Returns:
        F1 score per factor and over all factors.
    """
    if subset is None:
        subset = data["cell type"].unique()
    if flavor != "f1":
        assert correct_celltypes is not None, "correct_celltypes must be provided."

    if correct_celltypes is None:
        f1_stats = pd.DataFrame(columns=list(set(data.drop(columns=["cell type"]).columns)&set(subset)),
                                index=data["cell type"].unique()
                                )

        for factor in list(set(data.drop(columns=["cell type"]).columns)&set(subset)):
            for celltype in data["cell type"].unique():
                TP = data[factor][data["cell type"] == celltype].sum()
                FP = data[factor][data["cell type"] != celltype].sum()
                if general_stats is not None:
                    FN = general_stats.loc["factor_area_weighted" if weighted else "factor_area", factor] - TP - FP
                else:
                    FN = 1 - TP - FP
                f1_stats.loc[celltype, factor] = _f1_score(TP, FP, FN)
        f1_stats.sort_index(inplace=True)
        return f1_stats

    rates = {}

    for cluster in list(set(correct_celltypes.keys())&set(subset)):
        rates[cluster] = defaultdict(float)
        rates[cluster]["TP"] = data[cluster][[x in rates[cluster] for x in data["cell type"]]].sum()
        rates[cluster]["FP"] = data[cluster][[x not in rates[cluster] for x in data["cell type"]]].sum()
        if flavor == "micro" or flavor == "all" or general_stats is not None:
            assert general_stats is not None, "general_stats must be provided. Important for micro F1 score."
            rates[cluster]["FN"] = general_stats.loc["factor_area_weighted" if weighted else "factor_area", cluster] - rates[cluster]["TP"] - rates[cluster]["FP"]
        else:
            rates[cluster]["FN"] = 1-rates[cluster]["TP"]-rates[cluster]["FP"]
        rates[cluster]["F1_score"] = _f1_score(rates[cluster]["TP"], rates[cluster]["FP"], rates[cluster]["FN"])

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
    if flavor == "macro" or flavor == "all":
        names.append("macro F1_score")
        stats_total.append(mean(stats_total))
    elif flavor == "micro" or flavor == "all":
        names.append("micro F1_score")
        f1_stats = defaultdict(float)
        for stats in rates.values():
            f1_stats["TP"] += stats["TP"]
            f1_stats["FP"] += stats["FP"]
            f1_stats["FN"] += stats["FN"]
        f1_stats["F1_score"] = _f1_score(f1_stats["TP"], f1_stats["FP"], f1_stats["FN"])
        stats_total.append(f1_stats["F1_score"])
    stats_total = pd.DataFrame(stats_total, index=names, columns=["F1_statistics"]).T
    return stats_total

def f1_adapted(data: pd.DataFrame, min_value: Union[pd.DataFrame, float], correct_celltypes: Dict[str, str]=None) -> pd.DataFrame:
    """
    Computes the F1 score for mean and var information.
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
        assert set(data.drop(columns=["cell type"]).columns) <= set(min_value.index), \
            "some cell types present in data are not represented in min_value DataFrame."

    if correct_celltypes is None:
        f1_stats = pd.DataFrame(columns=data.drop(columns=["cell type"]).columns, index=data["cell type"].unique())

        for factor in data.drop(columns=["cell type"]).columns:
            for celltype in data["cell type"].unique():
                TP = (data.loc[data["cell type"] == celltype, factor] > float(min_value.loc[factor].values)).sum()
                FP = (data.loc[data["cell type"] != celltype, factor] > float(min_value.loc[factor].values)).sum()
                FN = len(data) - TP - FP
                f1_stats.loc[celltype, factor] = _f1_score(TP, FP, FN)
        f1_stats.sort_index(inplace=True)
        return f1_stats

    rates = {}

    for cluster, elements in correct_celltypes.items():
        rates[cluster] = defaultdict(float)
        rates[cluster]["TP"] = (data.loc[[x in elements for x in data["cell type"]], cluster] > min_value[cluster]).sum()
        rates[cluster]["FP"] = (data.loc[[x not in elements for x in data["cell type"]], cluster] > min_value[cluster]).sum()
        rates[cluster]["FN"] = len(data) - rates[cluster]["TP"] - rates[cluster]["FP"]
        rates[cluster]["F1_score"] = _f1_score(rates[cluster]["TP"], rates[cluster]["FP"], rates[cluster]["FN"])

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
