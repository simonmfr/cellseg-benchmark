import argparse

from cellseg_benchmark.metrics import (
    compute_metric_for_all_methods,
    compute_negative_marker_purity,
    get_negative_markers,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute negative marker purity for all methods in cohort."
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
    results_name = f"marker_gene_metrics/negative_marker_purity_{suffix}.csv"
    # load markers
    neg_marker_mask, ratio_celltype = get_negative_markers(
        args.cohort, vascular_subset=args.subset_vascular_celltypes
    )
    # compute negative marker purity
    compute_metric_for_all_methods(
        compute_negative_marker_purity,
        results_name=results_name,
        neg_marker_mask_sc=neg_marker_mask,
        ratio_celltype_sc=ratio_celltype,
        **vars(args),
    )
    # plot_MECR_score(args.cohort, suffix, show=False)
