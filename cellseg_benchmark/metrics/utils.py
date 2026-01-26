import os
from pathlib import Path

from cellseg_benchmark._constants import BASE_PATH


def compute_metric_for_all_methods(
    metric_func, cohort, base_path=BASE_PATH, methods=None, **kwargs
):
    """Generic function to iterate over all methods to compute metric.

    Args:
        metric_func: function to call for every method. Mandatory arguments cohort, method, base_path will be filled automatically.
        cohort: cohort to compute metric for
        base_path: base path to the data
        methods: list of methods to iterate over. If None, will iterate over all methods
        kwargs: further keyword arguments passed to metric_func
    """
    if methods is None:
        data_path = Path(base_path) / "analysis" / cohort
        methods = os.listdir(data_path)
    for i, method in enumerate(methods):
        print(f"{i + 1}/{len(methods)} ", end="")
        metric_func(cohort=cohort, method=method, base_path=base_path, **kwargs)
