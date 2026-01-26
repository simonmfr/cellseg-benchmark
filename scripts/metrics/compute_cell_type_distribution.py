import argparse

from cellseg_benchmark.metrics import (
    compute_cell_type_distribution,
    compute_metric_for_all_methods,
    plot_cell_type_distribution,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute cell type distributions for all methods in cohort."
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
        "--adata_name",
        default="adata_integrated",
        help="Name of the adata file (either adata_integrated or adata_vascular_subtypes)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing results"
    )

    args = parser.parse_args()

    results_name = f"cell_type_metrics/cell_type_distribution_{args.adata_name}_{args.celltype_name}.csv"
    compute_metric_for_all_methods(
        compute_cell_type_distribution, results_name=results_name, **vars(args)
    )
    plot_cell_type_distribution(args.cohort, f"{args.adata_name}_{args.celltype_name}")
