from collections import defaultdict
from statistics import mean
from typing import Dict, List, Literal, Union

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
    general_stats: pd.DataFrame | Dict[str, pd.DataFrame] = None,
    flavor: Literal["f1", "micro", "macro", "all"] = "f1",
    celltype_name: str | None = "celltype",
    correct_celltypes: Dict[str, List[str]] = None,
    weighted: bool = False,
    subset: List[str] = None,
) -> pd.DataFrame:
    """Computes the F1 score given area information.

    Args:
        data: Data with factors and celltype annotation.
        general_stats: If provided, data must be total for FN computation. Otherwise, data is relative
        flavor: "Micro" --> add micro f1, "Macro" --> add macro f1, "All" --> add micro and macro f1.
        correct_celltypes: specify correct celltype per factor.
        weighted: Data based on weighted or unweighted area.
        subset: Subset of ficture factors to compute F1 score for.

    Returns:
        F1 score per factor and over all factors.
    """
    if isinstance(general_stats, pd.DataFrame):
        assert isinstance(data, pd.DataFrame), "general_stats must be DataFrame, since data is DataFrame."
        data = {"tmp": data}
        general_stats = {"tmp": general_stats}
    elif isinstance(general_stats, dict):
        assert isinstance(data, dict), "general_stats must be a dict, since data is dict."
    if subset is None:
        subset = set.intersection(*[set(x.drop(columns=[celltype_name]).columns) for x in data.values()])
    if flavor != "f1":
        assert correct_celltypes is not None, "correct_celltypes must be provided."
    if flavor in ["micro", "all"]:
        assert general_stats is not None, "general_stats must be provided."

    if correct_celltypes is None:
        cols = list(set.intersection(*[set(x.drop(columns=[celltype_name]).columns) for x in data.values()]) & set(subset))
        cols.sort()
        celltypes_unique = set()
        for val in data.values():
            val[celltype_name] = val[celltype_name].cat.add_categories("Unknown")
            val[celltype_name].fillna("Unknown", inplace=True)
            celltypes_unique.update(val[celltype_name].unique())
        celltypes_unique = list(celltypes_unique)
        celltypes_unique.sort()
        f1_stats = pd.DataFrame(columns=cols, index=celltypes_unique)

        for factor in cols:
            for celltype in list(celltypes_unique):
                TP, FP, FN = 0, 0, 0
                for key, val in data.items():
                    TP_n = val.loc[val[celltype_name] == celltype, factor].sum()
                    FP_n = val.loc[val[celltype_name] != celltype, factor].sum()
                    if general_stats is not None:
                        FN += (
                            general_stats[key].loc[
                                "factor_area_weighted" if weighted else "factor_area",
                                factor,
                            ]
                            - TP_n
                            - FP_n
                        )
                    else:
                        FN += 1 - TP_n - FP_n
                    TP += TP_n
                    FP += FP_n
                f1_stats.loc[celltype, factor] = _f1_score(TP, FP, FN)
        f1_stats.sort_index(inplace=True)
        return f1_stats

    rates = {}
    cols = list(set(correct_celltypes.keys()) & set(subset))
    cols.sort()

    for cluster in cols:
        rates[cluster] = defaultdict(float)
        for key, val in data.items():
            TP_n = val.loc[
                [x in correct_celltypes[cluster] for x in val[celltype_name]], cluster
            ].sum()
            FP_n = val.loc[
                [x not in correct_celltypes[cluster] for x in val[celltype_name]], cluster
            ].sum()
            if flavor == "micro" or flavor == "all" or general_stats is not None:
                assert general_stats is not None, (
                    "general_stats must be provided. Important for micro F1 score."
                )
                rates[cluster]["FN"] = (
                    general_stats[key].loc[
                        "factor_area_weighted" if weighted else "factor_area", cluster
                    ]
                    - TP_n
                    - FP_n
                )
            else:
                rates[cluster]["FN"] = 1 - TP_n - FP_n
            rates[cluster]["TP"] += TP_n
            rates[cluster]["FP"] += FP_n
        rates[cluster]["F1_score"] = _f1_score(
            rates[cluster]["TP"], rates[cluster]["FP"], rates[cluster]["FN"]
        )

    stats_total = []
    names = []
    for cluster, stats in rates.items():
        names.append(cluster)
        stats_total.append(stats["F1_score"])
    if flavor == "macro" or flavor == "all":
        names.append("macro F1_score")
        stats_total.append(mean(stats_total))
    if flavor == "micro" or flavor == "all":
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
