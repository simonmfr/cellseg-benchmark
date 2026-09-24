#!/usr/bin/env python
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))

from cellseg_benchmark import BASE_PATH
from cellseg_benchmark.metrics import (
    compute_cell_density,
    compute_metric_for_all_methods,
    plot_cell_density,
    tissue_polygons,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute post-QC cells per mm² of tissue for all methods in cohort."
    )
    parser.add_argument("cohort", help="Cohort name.")
    parser.add_argument(
        "--methods",
        nargs="+",
        help="Methods to compute score for. If not specified, use all methods.",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing results"
    )
    args = parser.parse_args()
    tissue = tissue_polygons(args.cohort)
    compute_metric_for_all_methods(
        compute_cell_density, results_name="cell_density/cell_density.csv", tissue=tissue, **vars(args)
    )
    plot_cell_density(
        args.cohort, tissue, pathlib.Path(BASE_PATH) / "metrics" / args.cohort / "cell_density" / "cell_density.png"
    )
