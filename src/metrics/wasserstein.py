import pandas as pd
from scipy.stats import wasserstein_distance


def calculate_wasserstein_distance(
    adata, target_cluster, groupby="leiden", layer="librarysize_log1p_norm"
):
    """Calculate wasserstein distance."""
    cluster_mask = adata.obs[groupby] == target_cluster
    gene_expr = adata.layers[layer]

    wasserstein_distances = []
    for i, gene in enumerate(adata.var_names):
        expr_cluster = gene_expr[cluster_mask, i]
        expr_other = gene_expr[~cluster_mask, i]
        distance = wasserstein_distance(expr_cluster, expr_other)
        wasserstein_distances.append((gene, distance))

    distance_df = pd.DataFrame(wasserstein_distances, columns=["gene", "distance"])
    return distance_df.sort_values(by="distance", ascending=False).reset_index(
        drop=True
    )
