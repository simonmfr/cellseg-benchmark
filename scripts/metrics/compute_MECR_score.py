import argparse

from cellseg_benchmark.metrics import (
    compute_MECR_score,
    compute_metric_for_all_methods,
    plot_MECR_score,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute MECR score for all methods in cohort."
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
    parser.add_argument(
        "--subset_vascular_celltypes",
        action="store_true",
        help="Subset gene list to vascular celltypes",
    )

    args = parser.parse_args()
    suffix = "all"
    if args.subset_vascular_celltypes:
        suffix = "vascular_celltypes"
    results_name = f"marker_gene_metrics/MECR_score_{suffix}.csv"
    compute_metric_for_all_methods(
        compute_MECR_score, results_name=results_name, **vars(args)
    )
    plot_MECR_score(args.cohort, suffix, show=False)
