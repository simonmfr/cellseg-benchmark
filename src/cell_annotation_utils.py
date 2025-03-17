import os

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from matplotlib.pyplot import rc_context
from scipy.stats import median_abs_deviation
from spatialdata import SpatialData


def export_filter_adatas_from_sdata(sdata_path, sample_name, run_filter=True, min_counts=20, min_genes=5):
    """Extracts all adata objects from sdata.tables, filters cells based on given thresholds, and exports object into folder "_cell_type_annotation"."""
    sdata = SpatialData.read(sdata_path)
    
    # Create target dir
    export_dir = os.path.join(sdata_path, "..", "results", "_cell_type_annotation")
    os.makedirs(export_dir, exist_ok=True)
    
    for key, adata in sdata.tables.items():
        # Check if file already exists
        export_path = os.path.join(export_dir, f"{key}.h5ad")
        if os.path.exists(export_path):
            print(f"Skipping {key}: File already exists at {export_path}")
            continue
            
        print(f"Exporting {key}...")
        
        if not isinstance(adata, anndata.AnnData):
            print(f"Skipping {key}: Not an AnnData object.")
            continue
        
        # Cell filtering
        if run_filter:
            print(f"# cells before filtering: {adata.n_obs}")
            sc.pp.filter_cells(adata, min_counts=min_counts)
            sc.pp.filter_cells(adata, min_genes=min_genes)
            print(f"# cells after filtering: {adata.n_obs}")
        
        # Export adata
        adata.write_h5ad(export_path)


def flag_low_quality_mappings(df, group_col, value_col, mad_factor=5, inplace=True):
    """Flags low quality data by replacing group_col value with "Undefined" for rows with values significantly below the median (based on median absolute deviation).
    
    Args:
        df: DataFrame containing the data to analyze
        group_col: Column name used for grouping
        value_col: Column name containing values to evaluate
        mad_factor: Multiplier for median absolute deviation threshold (default: 5)
        inplace: If True, modifies df directly; if False, returns a copy (default: True)
       
    Returns:
        None if inplace=True, modified DataFrame copy if inplace=False
    """
    if inplace:
        low_quality_mask = df.groupby(group_col).apply(
        lambda x: x[value_col] < (x[value_col].median() - mad_factor * median_abs_deviation(x[value_col]))
        )
        df.loc[low_quality_mask.values, group_col] = "Undefined"
        return
    else:
        return df.copy()


def group_cell_types(meta):
    """Format and collapse cell types from MapMyCells into larger groups."""
    meta["cell_type"] = meta["allen_SUBC"].str.split().str[-1]

    neuron_mask = meta["cell_type"].isin(['Gaba', 'Glut', 'Gly-Gaba', 'Dopa-Gaba',
                                              'Hist-Gaba', 'Gaba-Glut', 'Dopa', 'Glut-Sero', 'Gaba-Chol', 'Chol',
                                              'Glut-Chol', 'Glyc-Gaba'])
    meta.loc[neuron_mask, "cell_type"] = "Neurons-" + meta.loc[neuron_mask, "cell_type"]

    nn_mask = meta["cell_type"] == "NN"
    meta.loc[nn_mask, "cell_type"] = meta.loc[nn_mask, "allen_SUBC"]

    meta["cell_type"] = meta["cell_type"].str.replace(r" NN$", "", regex=True)
    meta["cell_type"] = meta["cell_type"].str.replace(r"^\d+\s*", "", regex=True)

    replacement_dict = {
        "Astro-TE": "Astrocytes", "Astro-NT": "Astrocytes",
        "Astro-OLF": "Astrocytes", "Astro-CB": "Astrocytes", "Bergmann": "Astrocytes",
        "Oligo": "Oligodendrocytes",
        "Peri": "Pericytes", "Endo": "ECs", "VLMC": "VLMCs", "ABC": "VLMCs", "SMC": "SMCs",
        "OPC": "OPCs", "OEC": "OECs",
        "Astroependymal": "Ependymal", "CHOR": "Choroid Plexus", "Tanycyte": "Ependymal", "Hypendymal": "Ependymal",
        "Monocytes": "Immune-Other", "DC": "Immune-Other", "Lymphoid": "Immune-Other", "BAM": "BAMs",
        'IMN': 'Neurons-Immature',
        'Neurons-Glut-Chol': 'Neurons-Other',
        'Neurons-Gaba-Glut': 'Neurons-Other',
        'Neurons-Chol': 'Neurons-Other',
        'Neurons-Hist-Gaba': 'Neurons-Other',
        'Neurons-Gaba-Chol': 'Neurons-Gaba',
        'Neurons-Glut-Sero': 'Neurons-Other',
        'Neurons-Dopa-Gaba': 'Neurons-Other',
        'Neurons-Glyc-Gaba': 'Neurons-Other',
        'Neurons-Gly-Gaba': 'Neurons-Other'
    }
    meta["cell_type"] = meta["cell_type"].replace(replacement_dict)

    return meta


def normalize_counts(adata, method="area", target_sum=250):
    """Normalize counts by area/volume or library size, then compute log1p and z-score.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with raw counts in .layers["counts"]
    method : str, default="area"
        Normalization method: "area", "volume", or "library"
    target_sum : int, default=250
        Target sum of counts per cell after normalization (250 used by Allen et al., 2024, Cell)
        
    Returns:
    --------
    None, modifies adata in-place
    """
    if method == "library":
        # Library size normalization
        adata.layers["library_log1p_norm"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=target_sum, layer="library_log1p_norm")
        sc.pp.log1p(adata, layer="library_log1p_norm")
        return adata
    
    # Area/volume normalization
    if method == "area":
        if "area" in adata.obs.columns:
            size_factor = adata.obs["area"].values
            norm_layer_name = "area_norm"
        elif "volume" in adata.obs.columns:
            print("Area not found, using volume instead.")
            size_factor = adata.obs["volume"].values
            norm_layer_name = "volume_norm"
        else:
            raise ValueError("Neither 'area' nor 'volume' found in adata.obs.columns")
    elif method == "volume":
        if "volume" in adata.obs.columns:
            size_factor = adata.obs["volume"].values
            norm_layer_name = "volume_norm"
        else:
            raise ValueError("'volume' not found in adata.obs.columns")
    else:
        raise ValueError("Method must be one of: 'area', 'volume', or 'library'")
        
    # Normalize by area/volume
    adata.layers[norm_layer_name] = sp.csr_matrix(adata.X, dtype=np.float32).multiply(1 / size_factor[:, None])
    
    # Scale to target sum
    norm_matrix = adata.layers[norm_layer_name]
    scaling_factors = target_sum / norm_matrix.sum(axis=1).A1
    adata.layers[norm_layer_name] = norm_matrix.multiply(scaling_factors[:, None])
    
    # Verify mean count matches target
    assert round(adata.layers[norm_layer_name].sum(axis=1).A1.mean()) == target_sum
    
    # Log-transform and z-score
    log_layer_name = f"{norm_layer_name.split('_')[0]}_log1p_norm"
    adata.layers[log_layer_name] = sc.pp.log1p(adata.layers[norm_layer_name], copy=True)
    adata.layers["zscore"] = sc.pp.scale(adata.layers[log_layer_name], zero_center=True, max_value=None, copy=True)
    
    return adata


def plot_mad_thresholds(Allen_MMC_metadata, out_path, name = "mad_threshold", group_column="allen_CLUS", value_column="allen_avg_cor_CLUS", mad_factor=3, figsize=(13, 7)):
    """Plot a violin plot of data grouped by a specified column, with horizontal lines showing MAD thresholds.

    Parameters:
    - Allen_MMC_metadata: DataFrame containing the data
    - group_column: Column name for grouping (e.g. "allen_CLUS")
    - value_column: Column name for the values (e.g. "allen_avg_cor_CLUS")
    - mad_factor: The factor for excluding outliers based on MAD (default is 3)
    """
    # Compute MAD within each cluster
    grouped_mad = Allen_MMC_metadata.groupby(group_column)[value_column].apply(median_abs_deviation)
    unique_groups = Allen_MMC_metadata[group_column].unique()

    with rc_context({'figure.figsize': figsize}):
        plt.grid(True, zorder=0)

        # Create violin plot
        sns.violinplot(data=Allen_MMC_metadata, x=group_column, y=value_column,
                       density_norm="width", inner="quart", zorder=2)

        # Add horizontal lines for each group's MAD threshold
        for group, mad in grouped_mad.items():
            group_median = Allen_MMC_metadata[Allen_MMC_metadata[group_column] == group][value_column].median()
            threshold = group_median - mad_factor * mad

            plt.axhline(y=threshold, color='r', linestyle='-',
                        xmin=(list(unique_groups).index(group) - 0) / len(unique_groups),
                        xmax=(list(unique_groups).index(group) + 0.9) / len(unique_groups),
                        zorder=2)

        plt.xticks(rotation=90, ha='center')
        plt.tight_layout()
        plt.savefig(os.path.join(out_path, name + '.png'))
        #plt.close()
        
        return

    
def plot_mapping_qc(json_results, adata, mapping_result, level='CCN20230722_CLAS'):
    """Create quality control plots for cell type mapping results.

    Parameters:
    -----------
    json_results : dict
        Results from cell type mapping, containing taxonomy tree and mapping results
    adata : anndata.AnnData
        Annotated data matrix with cell information
    level : str, optional
        Hierarchical level to analyze (default is 'CCN20230722_CLAS')

    Returns:
    --------
    fig : matplotlib.figure.Figure
        Figure containing three subplots of mapping QC metrics
    """
    # Get the level name from taxonomy tree
    level_name = json_results['taxonomy_tree']['hierarchy_mapper'][level]

    # Get cell barcodes and create mapping result dictionary
    cell_barcodes = adata.obs.index.values

    # Calculate number of non-zero genes for each cell
    n_genes = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, "A1") else (adata.X > 0).sum(axis=1)

    # Extract average correlation and bootstrapping probability
    avg_corr = np.array([mapping_result[b][level]['avg_correlation'] for b in cell_barcodes])
    bootstrapping_prob = np.array([mapping_result[b][level]['bootstrapping_probability'] for b in cell_barcodes])

    # Create figure with three subplots
    fig = plt.figure(figsize=(15, 4))

    # Subplot 1: Number of genes vs Average Correlation
    corr_axis = fig.add_subplot(1, 3, 1)
    corr_axis.set_title(level_name)
    corr_axis.scatter(n_genes, avg_corr, s=1)
    corr_axis.set_xlabel('N non-zero genes')
    corr_axis.set_ylabel('Average Correlation')
    corr_axis.set_xscale('log')

    # Subplot 2: Number of genes vs Bootstrapping Probability
    prob_axis = fig.add_subplot(1, 3, 2)
    prob_axis.scatter(n_genes, bootstrapping_prob, s=1)
    prob_axis.set_xlabel('N non-zero genes')
    prob_axis.set_ylabel('Bootstrapping Probability')
    prob_axis.set_xscale('log')

    # Subplot 3: Average Correlation vs Bootstrapping Probability
    corr_v_prob_axis = fig.add_subplot(1, 3, 3)
    corr_v_prob_axis.scatter(avg_corr, bootstrapping_prob, s=1)
    corr_v_prob_axis.set_xlabel('Average Correlation')
    corr_v_prob_axis.set_ylabel('Bootstrapping Probability')

    return fig


def process_allen_metadata(json_results, mapping_result):
    """Process Allen metadata from JSON results and mapping results.

    Parameters:
    -----------
    json_results : dict
        Dictionary containing taxonomy tree information
    mapping_result : dict
        Dictionary containing mapping results for each barcode

    Returns:
    --------
    pd.DataFrame
        Processed metadata DataFrame with readable labels and correlations
    """
    # Build a mapping of node to human-readable names for all levels
    node_to_label = {
        level: {node: json_results['taxonomy_tree']['name_mapper'][level][node]['name']
                for node in json_results['taxonomy_tree'][level]}
        for level in json_results['taxonomy_tree']['hierarchy']
    }

    data = {}
    for level in json_results['taxonomy_tree']['hierarchy']:
        for barcode, result in mapping_result.items():
            if barcode not in data:
                data[barcode] = {}

            if level in result:
                assignment = result[level]['assignment']
                level_short = level.replace('CCN20230722_', '')  # Remove 'CCN20230722' from level name

                data[barcode][f'label_id_{level_short}'] = assignment
                data[barcode][f'allen_avg_cor_{level_short}'] = result[level]['avg_correlation']
                data[barcode][f'allen_{level_short}'] = node_to_label[level].get(assignment, 'Unknown')

    # Convert to DataFrame
    Allen_MMC_metadata = pd.DataFrame.from_dict(data, orient='index')

    # Drop columns with non-readable names (label_id columns)
    Allen_MMC_metadata.drop(columns=Allen_MMC_metadata.filter(like="label_id_"), inplace=True)

    return Allen_MMC_metadata


# Define matching dictionary for cell types (MapMyCells/Yao2023 cell annotation) to score columns (scoring of Yao2023 marker genes)
match_dict = {
   "Astrocytes": "Astrocytes",
   "BAMs": "BAMs",
   "Choroid Plexus": "Ependymal",
   "ECs": "ECs",
   "Ependymal": "Ependymal", 
   "Immune-Other": "Immune-Other",
   "Microglia": "Microglia",
   "Neurons-Dopa": "Neurons-Dopa",
   "Neurons-Glut": "Neurons-Glut",
   "Neurons-Other": "Neurons-Other",
   "Neurons-Immature": "Neurons-Immature",
   "Neurons-Gaba": "Neurons-Gaba",
   "OECs": "OECs",
   "OPCs": "OPCs",
   "Oligodendrocytes": "Oligodendrocytes",
   "Pericytes": "Pericytes",
   "SMCs": "SMCs",
   "VLMCs": "VLMCs"
}

# 1. Crosstab calculation for cell type assignment
def assign_cell_types_by_cluster(adata, leiden_res=3.0, min_cells=100):
    """Assign cell types to leiden clusters based on majority vote."""
    leiden_col = f"leiden_res{leiden_res}".replace(".", "_")
    cell_type_col = "cell_type"
    
    # Create crosstab
    ctab = pd.crosstab(
        adata.obs[leiden_col],
        adata.obs[cell_type_col],
        margins=True,
        margins_name="Total"
    )
    
    # Exclude cell types with < min_cells cells
    ctab = ctab.loc[:, ctab.loc["Total"] >= min_cells]
    
    # Normalize column-wise
    ctab_percentage = ctab.div(ctab.loc["Total"], axis=1) * 100
    
    # Get cluster labels based on max percentage
    cluster_labels = ctab_percentage.drop(index="Total", columns="Total").idxmax(axis=1)
    
    # Convert to dictionary
    cluster_labels_dict = cluster_labels.to_dict()
    
    # Assign cell types to clusters
    adata.obs["cell_type_clusters"] = adata.obs[leiden_col].map(cluster_labels_dict)
    
    return adata, cluster_labels_dict, ctab_percentage

# 2. Score cell types using marker genes
def score_cell_types(adata, marker_genes_dict, top_n_genes=25):
    """Score cells for marker gene expression using only top n genes.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    marker_genes_dict : dict
        Dictionary of cell types and their marker genes
    top_n_genes : int, default=30
        Number of top genes to use for scoring
    """
    for cell_type, genes in marker_genes_dict.items():
        if genes:  # Only score if gene list is not empty
            # Take only the top n genes
            genes_to_use = genes[:top_n_genes] if len(genes) > top_n_genes else genes
            score_name = f"score_{cell_type}"
            sc.tl.score_genes(adata, genes_to_use, score_name=score_name)
    
    return adata

# 3. Final cell type assignment based on scores
def assign_final_cell_types(adata, cluster_labels_dict, match_dict, leiden_res=3.0, score_threshold=0.25):
    """Assign final cell types based on scoring results."""
    leiden_col = f"leiden_res{leiden_res}".replace(".", "_")
    
    # Initialize final cell types with Undefined
    adata.obs["cell_type_final"] = "Undefined"
    
    # Store reasons for undefined assignments
    undefined_reasons = {}
    
    # Initialize cluster-score matrix
    clusters = list(adata.obs[leiden_col].unique())
    cell_types = list(match_dict.keys())
    cluster_score_matrix = pd.DataFrame(index=clusters, columns=cell_types, dtype=float)
    
    # For each cluster
    for cluster in clusters:
        # Get the cell type assigned by crosstab
        assigned_cell_type = cluster_labels_dict.get(cluster, "Unknown")
        
        # Get cells in this cluster
        mask = adata.obs[leiden_col] == cluster
        
        # Collect all scores for this cluster
        all_scores = {}
        for cell_type, score_suffix in match_dict.items():
            score_col = f"score_{score_suffix}"
            if score_col in adata.obs.columns:
                mean_score = adata.obs.loc[mask, score_col].mean()
                all_scores[cell_type] = mean_score
                cluster_score_matrix.loc[cluster, cell_type] = mean_score
            else:
                cluster_score_matrix.loc[cluster, cell_type] = np.nan
        
        if not all_scores:
            undefined_reasons[cluster] = "No scores found for any cell type"
            continue
            
        # Find highest scoring cell type
        highest_cell_type = max(all_scores, key=all_scores.get)
        highest_score = all_scores[highest_cell_type]
        
        # Get score for assigned cell type
        assigned_score_col = f"score_{match_dict.get(assigned_cell_type, assigned_cell_type)}"
        if assigned_score_col in adata.obs.columns:
            assigned_score = adata.obs.loc[mask, assigned_score_col].mean()
        else:
            undefined_reasons[cluster] = f"No score column found for {assigned_cell_type}"
            continue
        
        # If highest score is above threshold, use it
        if highest_score > score_threshold:
            if highest_cell_type != assigned_cell_type:
                print(f"Cluster {cluster}: Reassigned from {assigned_cell_type} to {highest_cell_type} (higher score: {highest_score:.4f} vs {assigned_score:.4f})")
            adata.obs.loc[mask, "cell_type_final"] = highest_cell_type
        else:
            undefined_reasons[cluster] = f"No cell type has score above threshold {score_threshold} (highest: {highest_cell_type} with {highest_score:.4f})"
    
    return adata, undefined_reasons, cluster_score_matrix

# Main function to run the entire pipeline
def cell_type_annotation_pipeline(adata, marker_genes_dict, leiden_res=3.0, min_cells=100, score_threshold=0.25, top_n_genes=25):
    """Run the complete cell type annotation pipeline."""
    # Check if leiden clustering exists, run if not present
    leiden_col = f"leiden_res{leiden_res}".replace(".", "_")
    if leiden_col not in adata.obs.columns:
        print(f"Running leiden clustering with resolution {leiden_res}...")
        sc.tl.leiden(adata, key_added=leiden_col, resolution=leiden_res)
    
    print("1. Assigning cell types based on cluster composition...")
    adata, cluster_labels_dict, normalized_crosstab = assign_cell_types_by_cluster(adata, leiden_res, min_cells)
    
    print("2. Scoring cells for marker gene expression...")
    adata = score_cell_types(adata, marker_genes_dict, top_n_genes)
    
    print("3. Assigning final cell types based on scores...")
    adata, undefined_reasons, cluster_score_matrix = assign_final_cell_types(adata, cluster_labels_dict, match_dict, leiden_res, score_threshold)
    
    # Print undefined reasons
    print("\nClusters marked as Undefined:")
    for cluster, reason in undefined_reasons.items():
        print(f"Cluster {cluster}: {reason}")
    
    # Calculate summary statistics
    defined_count = (adata.obs["cell_type_final"] != "Undefined").sum()
    total_count = len(adata.obs)
    defined_percent = defined_count / total_count * 100
    
    print(f"\nAnnotation summary: {defined_count}/{total_count} cells ({defined_percent:.2f}%) assigned to defined cell types")
    
    # Create summary dataframe
    summary_df = pd.DataFrame({
        'Cluster': list(cluster_labels_dict.keys()),
        'Initial_Cell_Type': [cluster_labels_dict[cluster] for cluster in cluster_labels_dict.keys()],
        'Final_Cell_Type': [adata.obs[adata.obs[f"leiden_res{leiden_res}".replace('.', '_')] == cluster]['cell_type_final'].iloc[0] 
                           if len(adata.obs[adata.obs[f"leiden_res{leiden_res}".replace('.', '_')] == cluster]) > 0 else 'Undefined' 
                           for cluster in cluster_labels_dict.keys()]
    })
    
    return adata, summary_df, normalized_crosstab, cluster_score_matrix