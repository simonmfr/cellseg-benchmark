import argparse

from cellseg_benchmark.metrics import (
    compute_metric_for_all_methods,
    extract_general_stats,
    plot_general_stats,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute general stats for all methods in cohort."
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
    results_name = "general_stats/general_stats.csv"

    # compute negative marker purity
    compute_metric_for_all_methods(
        extract_general_stats,
        results_name=results_name,
        **vars(args),
    )
    for metric in [
        "volume_final",
        "area",
        "sphericity",
        "elongation",
        "intensities_PolyT",
        "intensities_DAPI",
        "Ovrlpy_stats_mean_integrity",
    ]:
        plot_general_stats(args.cohors, metric, show=False)
