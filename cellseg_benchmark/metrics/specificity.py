import numpy as np
import scanpy as sc
import scipy.sparse as sp


def scanpy_volume_log_norm(
    adata,
    volume_key,
    target_sum=250,
    layer_norm_name="volume_norm",
    layer_log_name="volume_log1p_norm",
):
    """Performs volume normalization, scaling, and log transformation on AnnData object.

    Args:
        adata: AnnData object containing count data in .X
        volume_key: Key in adata.obs containing cell volume measurements
        target_sum: Target sum for scaling each cell (default 250 as in Allen et al. 2024)
        layer_norm_name: Name for normalized layer (default "volume_norm")
        layer_log_name: Name for log-transformed layer (default "volume_log1p_norm")

    Returns:
        Modified AnnData object with new layers
    """
    if not np.allclose(adata.X.max(), round(adata.X.max())):
        raise ValueError("adata.X does not appear to be count data")

    # Step 1: Normalize by cell volume
    volumes = adata.obs[volume_key].values
    sparse_matrix = sp.csr_matrix(adata.X, dtype=np.float32)
    normalized_matrix = sparse_matrix.multiply(1 / volumes[:, None])
    adata.layers[layer_norm_name] = normalized_matrix

    # Step 2: Scale each cell's total expression to sum to target_sum
    row_sums = normalized_matrix.sum(axis=1).A1
    scaling_factors = target_sum / row_sums
    normalized_matrix = normalized_matrix.multiply(scaling_factors[:, None])
    normalized_matrix = sp.csr_matrix(normalized_matrix, dtype=np.float32)
    adata.layers[layer_norm_name] = normalized_matrix

    # Verify
    row_sums = normalized_matrix.sum(axis=1).A1
    assert round(row_sums.mean()) == 250

    # Step 3: Log-transform
    adata.layers[layer_log_name] = sc.pp.log1p(adata.layers[layer_norm_name], copy=True)

    return adata


def scanpy_pca_umap_leiden(
    adata,
    leiden_res=[0.5, 1.0],
    layer="volume_log1p_norm",
    hvg=False,
    hvg_n_top_genes=None,
):
    """Perform PCA, UMAP, and Leiden clustering on an AnnData object.

    Parameters:
    - adata: AnnData
        Anndata object.
    - leiden_res: list, default=[0.5, 1.0]
        List of resolution values for Leiden clustering.
    - layer: str, default="volume_log1p_norm"
        The layer to use as input for analysis.
    - hvg: bool, default=False
        Whether to select highly variable genes before PCA.
    - hvg_n_top_genes: int, optional
        Number of top highly variable genes to select if hvg=True.
    """
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found in adata.layers")
    adata.X = adata.layers[layer]
    if hvg:
        sc.pp.highly_variable_genes(adata, n_top_genes=hvg_n_top_genes)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    for res in leiden_res:
        sc.tl.leiden(
            adata, key_added=f"leiden_res{res}".replace(".", "_"), resolution=res
        )
    return adata


def calculate_marker_specificity_sensitivity(
    adata, target_marker, groupby, target_cluster, layer="counts", zero_threshold=0
):
    """Quantifying False Positive Rate/Specificity/etc. by calculating the ratio of the expression level of the marker in the cluster of interest to its expression in other clusters."""
    adata.X = adata.layers[layer]
    target_marker_adata = adata[:, adata.var_names.str.match(target_marker)]

    cells_cluster = adata[adata.obs[groupby] == target_cluster].obs_names
    cells_other = adata[adata.obs[groupby] != target_cluster].obs_names
    adata_cluster = target_marker_adata[cells_cluster]
    adata_other = target_marker_adata[cells_other]

    # Calculate metrics
    TP = (adata_cluster.X > zero_threshold).sum()  # cluster expressing
    FN = (adata_cluster.X <= zero_threshold).sum()  # cluster nonexpressing
    FP = (adata_other.X > zero_threshold).sum()  # other expressing
    TN = (adata_other.X <= zero_threshold).sum()  # other nonexpressing

    metrics = {
        "marker": target_marker,
        "group": groupby,
        "target_cluster": target_cluster,
        "expression_threshold": zero_threshold,
        "metrics": {
            "FPR": round(FP / (FP + TN), 3),
            "specificity": round(TN / (TN + FP), 3),
            "sensitivity": round(TP / (TP + FN), 3),
            "PPV": round(TP / (TP + FP), 3),
            "NPV": round(TN / (TN + FN), 3),
        },
    }

    return metrics
