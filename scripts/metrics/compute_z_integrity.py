#!/usr/bin/env python
"""Compute the z integrity metric for a cohort's 2D segmentation methods and save to CSV.

Same as the compute cell of notebooks/metrics_z_integrity.ipynb. Needs the ovrlpy z integrity
maps (scripts/seg_postprocessing/ovrlpy_z_integrity_map.py) and the integrated adatas.
"""

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))

from cellseg_benchmark import BASE_PATH, _constants
from cellseg_benchmark.metrics import (
    compute_metric_for_all_methods,
    compute_z_integrity,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohort", help="Cohort name.")
    parser.add_argument("--methods", nargs="+", help="Methods to score (default: all methods in the cohort).")
    parser.add_argument("--subset_vascular_celltypes", action="store_true", help="Score only the cells of the vascular subset.")
    parser.add_argument("--overwrite", action="store_true", help="Recompute methods already present in the CSV.")
    args = parser.parse_args()

    vascular = args.subset_vascular_celltypes
    methods = args.methods or sorted(p.name for p in (pathlib.Path(BASE_PATH) / "analysis" / args.cohort).iterdir())
    methods = [m for m in methods if not m.startswith(tuple(_constants.methods_3D))]
    compute_metric_for_all_methods(
        compute_z_integrity,
        cohort=args.cohort,
        results_name=f"z_integrity/z_integrity_{'vascular_celltypes' if vascular else 'all'}.csv",
        adata_name="adata_vascular_subset" if vascular else "adata_integrated",
        methods=methods,
        overwrite=args.overwrite,
    )
