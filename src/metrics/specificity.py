import scanpy as sc
import numpy as np

def calculate_marker_specificity_sensitivity(adata, target_marker, groupby, target_cluster, layer="counts", zero_threshold=0):
    """
    Quantifying False Positive Rate/Specificity/etc. by calculating the ratio of the expression level of the marker in the cluster of interest to its expression in other clusters.

    """
    adata.X = adata.layers[layer]
    target_marker_adata = adata[:, adata.var_names.str.match(target_marker)]
    
    cells_cluster = adata[adata.obs[groupby] == target_cluster].obs_names
    cells_other = adata[adata.obs[groupby] != target_cluster].obs_names
    adata_cluster = target_marker_adata[cells_cluster]
    adata_other = target_marker_adata[cells_other]
    
    # Calculate metrics
    TP = (adata_cluster.X > zero_threshold).sum() # cluster expressing
    FN = (adata_cluster.X <= zero_threshold).sum() # cluster nonexpressing
    FP = (adata_other.X > zero_threshold).sum() # other expressing
    TN = (adata_other.X <= zero_threshold).sum() # other nonexpressing
    
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
