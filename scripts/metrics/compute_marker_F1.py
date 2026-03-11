import argparse

from cellseg_benchmark.metrics import (
    compute_marker_F1_score,
    compute_metric_for_all_methods,
    plot_marker_F1_score,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute marker F1 score for all methods in cohort."
    )
    parser.add_argument("cohort", help="Cohort name.")
    parser.add_argument(
        "--methods",
        nargs="+",
        help="Methods to compute score for. If not specified, use all methods.",
    )
    parser.add_argument(
        "--celltype_name",
        default="cell_type_revised",
        help="Nme of celltype column. Default is cell_type_revised.",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing results"
    )

    args = parser.parse_args()
    results_name = f"marker_gene_metrics/marker_f1_score_{args.celltype_name}.csv"
    compute_metric_for_all_methods(
        compute_marker_F1_score, results_name=results_name, **vars(args)
    )

    plot_marker_F1_score(args.cohort, args.celltype_name, show=False)
