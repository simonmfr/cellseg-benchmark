import argparse

from cellseg_benchmark.metrics import (
    compute_metric_for_all_methods,
    compute_assigned_transcripts,
    plot_assigned_transcripts,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute assigned transcripts for all methods in cohort, splitting data per sample."
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
    results_name = "assigned_transcripts/assigned_transcript_counts.csv.csv"
    compute_metric_for_all_methods(
        compute_assigned_transcripts, results_name=results_name, pass_method=True, **vars(args)
    )
    plot_assigned_transcripts(args.cohort)