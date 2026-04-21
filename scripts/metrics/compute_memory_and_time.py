import argparse

from  cellseg_benchmark.metrics import (
    compute_metric_for_all_methods,
    extract_mem_and_time,
    plot_mem_and_time,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute memory and time usage for all methods in cohort."
    )
    parser.add_argument("cohort", help="Cohort name.")
    parser.add_argument(
        "plot_metric", choices=["memory", "cpus", "duration"], help="Specify the metric to plot."
    )
    parser.add_argument(
        "--ref_file_path", help="Path to file containing the raw output of sacct."
    )
    parser.add_argument(
        "--metrics_dir", help="Path to directory containing from sacct output extracted metrics."
    )
    args = parser.add_argument(
        "--ignore_missing", action="store_true", help="Ignore missing method metrics."
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing results"
    )

    args = parser.parse_args()
    results_name = f"Mem_and_time/mem_and_time.csv"

    #Ensure that the default values of ref_file_path and metrics_dir are used if arguments are not provided
    if args.ref_file_path is None:
        delattr(args, 'ref_file_path')
    if args.metrics_dir is None:
        delattr(args, 'metrics_dir')

    compute_metric_for_all_methods(
        extract_mem_and_time, results_name=results_name, **vars(args)
    )

    plot_mem_and_time(args.cohort, args.plot_metric, show=False)