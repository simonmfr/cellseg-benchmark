#!/usr/bin/env python
"""Compute the FICTURE F1 metric for a cohort's segmentation methods and save to CSV.

Results are written to ``metrics/<cohort>/ficture/ficture_f1.csv`` with one row per
(method, sample, cell type) plus pooled ``all_samples`` rows.
"""

import argparse

from cellseg_benchmark.metrics import compute_ficture_f1
from cellseg_benchmark.metrics.utils import compute_metric_for_all_methods

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("cohort", help="Cohort name.")
parser.add_argument(
    "--methods",
    nargs="+",
    help="Methods to score (default: all methods in the cohort).",
)
parser.add_argument(
    "--n_jobs", type=int, default=-1, help="Parallel workers over samples."
)
parser.add_argument(
    "--overwrite",
    action="store_true",
    help="Recompute methods already present in the CSV.",
)
args = parser.parse_args()

compute_metric_for_all_methods(
    compute_ficture_f1,
    cohort=args.cohort,
    results_name="ficture/ficture_f1.csv",
    methods=args.methods,
    n_jobs=args.n_jobs,
    overwrite=args.overwrite,
)
