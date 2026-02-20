import argparse

from cellseg_benchmark.metrics import (
    compute_clustering_scores,
    compute_metric_for_all_methods,
    plot_clustering_scores,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute clustering score (CH/SH) for all methods in cohort, splitting data per sample."
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
        help="Name of celltype column. Default is cell_type_revised. If set to 'leiden', leiden clustering is computed on the adata",
    )
    parser.add_argument(
        "--adata_name",
        default="adata_integrated",
        help="Name of the adata file (either adata_integrated or adata_vascular_subtypes)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing results"
    )
    parser.add_argument(
        "--sample_size",
        default=10000,
        help="Downsamples data to sample_size to speed up computation. Default is 10000.",
    )
    parser.add_argument(
        "--n_pcs",
        default=30,
        help="Number of principal components to use. Default is 30.",
    )
    parser.add_argument(
        "--leiden_resolution",
        default=1,
        help="Resolution to use for leiden clustering if desired. Default is 1",
    )

    args = parser.parse_args()
    results_name = (
        f"cell_type_metrics/clustering_score_{args.adata_name}_{args.celltype_name}.csv"
    )
    compute_metric_for_all_methods(
        compute_clustering_scores, results_name=results_name, **vars(args)
    )
    plot_clustering_scores(args.cohort, f"{args.adata_name}_{args.celltype_name}")
