import os
import sys
import subprocess
from datetime import date

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from matplotlib.pyplot import rc_context
from scipy.stats import median_abs_deviation


def assign_cell_types_to_clusters(adata, leiden_col, cell_type_col = "cell_type_mapmycells", min_cells=100):
    """Assign cell type labels to leiden clusters based on majority vote of MapMyCells results."""
    # Create crosstab
    cluster_cell_type_crosstab = pd.crosstab(
        adata.obs[leiden_col],
        adata.obs[cell_type_col],
        margins=True,
        margins_name="Total",
    )

    # Exclude cell types with fewer than the minimum number of cells
    cluster_cell_type_crosstab = cluster_cell_type_crosstab.loc[:, cluster_cell_type_crosstab.loc["Total"] >= min_cells]

    # Normalize column-wise to get percentages
    normalized_percentage = cluster_cell_type_crosstab.div(cluster_cell_type_crosstab.loc["Total"], axis=1) * 100

    # Identify the most frequent cell type for each cluster (maximum percentage)
    assigned_cell_type_dict = normalized_percentage.drop(index="Total", columns="Total").idxmax(axis=1).to_dict()

    return assigned_cell_type_dict, normalized_percentage


def assign_final_cell_types(
    adata, cluster_labels_dict, match_dict, leiden_col, out_col = "cell_type_final", score_threshold=0.25, logger=None
):
    """Assign final cell types based on scoring results in-place.
    
    Parameters:
    -----------
    adata : anndata.AnnData
        Anndata object
    cluster_labels_dict : dict
        Dictionary mapping cluster IDs to initially assigned cell types, created by assign_cell_types_by_cluster()
    match_dict : dict
        Dictionary mapping MapMyCells cell type names to their corresponding scanpy.tl.score_genes suffixes
    leiden_col : string
        Column in adata.obs containing Leiden clustering
    score_threshold : float, default=0.25
        Minimum score required to assign a cell type
    logger : logger, default=None
    
    Returns:
    --------
    adata : anndata.AnnData
        Updated anndata object with out_col in adata.obs
    undefined_reasons : dict
        Dictionary explaining why some clusters remained undefined
    cluster_score_matrix : pandas.DataFrame
        Matrix of scores for each cluster-cell type combination
    """

    # Initialize final cell types with Undefined
    adata.obs[out_col] = np.nan

    # Store reasons for undefined assignments
    undefined_reasons = {}

    # Initialize cluster-score matrix
    clusters = list(adata.obs[leiden_col].unique())
    cell_types = list(match_dict.keys())
    cluster_score_matrix = pd.DataFrame(index=clusters, columns=cell_types, dtype=float)

    # For each cluster
    for cluster in clusters:
        # Get the cell type assigned by crosstab
        assigned_cell_type = cluster_labels_dict.get(cluster, "Undefined")

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
        assigned_score_col = (
            f"score_{match_dict.get(assigned_cell_type, assigned_cell_type)}"
        )
        if assigned_score_col in adata.obs.columns:
            assigned_score = adata.obs.loc[mask, assigned_score_col].mean()
        else: # e.g. "Undefined" or "Mixed"
            assigned_score = np.nan

        # If highest score is above threshold, use it
        if highest_score > score_threshold:
            if highest_cell_type != assigned_cell_type:
                logger.warning(f"    Cluster {cluster}: Reassigned from {assigned_cell_type} to {highest_cell_type} (higher score: {assigned_score:.4f} vs {highest_score:.4f})")
            adata.obs.loc[mask, out_col] = highest_cell_type
        else:
            undefined_reasons[cluster] = (
                f"No cell type has score above threshold {score_threshold} (highest: {highest_cell_type} with {highest_score:.4f}). Set to undefined."
            )
            adata.obs.loc[mask, out_col] = "Undefined"

    return adata, undefined_reasons, cluster_score_matrix


def create_mixed_cell_types(df, 
                            diff_threshold=0.5, 
                            main_col='allen_SUBC', 
                            diff_col='allen_diff_rup1_prob_SUBC', 
                            runnerup_col='allen_runner_up_1_SUBC'
                           ):
   """Create DataFrame with mixed cell type annotations, based on MapMyCells probability differences."""
   # Identify mixed cells
   is_mixed = (df[diff_col] < diff_threshold) & (df[main_col] != df[runnerup_col])
   
   # Create result DataFrame
   result = pd.DataFrame(index=df.index)
   result['allen_SUBC_is_mixed'] = is_mixed.map({True: "mixed", False: "unique"})
   result['allen_SUBC_incl_mixed'] = df[main_col].where(~is_mixed, "Mixed")
   result['allen_SUBC_mixed_names'] = df[main_col].where(~is_mixed, df[main_col] + "/" + df[runnerup_col])
   
   return result
    

def export_filter_adatas_from_sdata(sdata_path, logger):
    """Extracts all adata objects from sdata.tables, filters cells based on given thresholds, and exports object into folder "_cell_type_annotation"."""

    from spatialdata import read_zarr
    
    sdata = read_zarr(sdata_path)

    # Create target dir
    export_dir = os.path.join(sdata_path, "..", "_cell_type_annotation")
    os.makedirs(export_dir, exist_ok=True)

    for key, adata in sdata.tables.items():
        # Check if file already exists
        export_path = os.path.join(export_dir, f"{key}.h5ad")
        if os.path.exists(export_path):
            logger.error(f"Skipping {key}: File already exists at {export_path}")
            continue

        logger.info(f"Extracting {key}")

        if not isinstance(adata, anndata.AnnData):
            logger.warning(f"Skipping {key}: AnnData object is not an AnnData object.")
            continue

        # Export adata
        adata.write_h5ad(export_path)


def group_cell_types(metadata_col):
    """Format and bin cell types from MapMyCells into defined groups."""

    result = metadata_col.str.split().str[-1].copy()

    neuron_mask = result.isin(
        [
            "Gaba",
            "Glut",
            "Gly-Gaba",
            "Dopa-Gaba",
            "Hist-Gaba",
            "Gaba-Glut",
            "Dopa",
            "Glut-Sero",
            "Gaba-Chol",
            "Chol",
            "Glut-Chol",
            "Glyc-Gaba",
        ]
    )
    result.loc[neuron_mask] = "Neurons-" + result.loc[neuron_mask]

    nn_mask = result == "NN"
    result.loc[nn_mask] = metadata_col.loc[nn_mask]

    result = result.str.replace(r" NN$", "", regex=True)
    result = result.str.replace(r"^\d+\s*", "", regex=True)

    replacement_dict = {
        "Astro-TE": "Astrocytes",
        "Astro-NT": "Astrocytes",
        "Astro-OLF": "Astrocytes",
        "Astro-CB": "Astrocytes",
        "Oligo": "Oligodendrocytes",
        "Peri": "Pericytes",
        "Endo": "ECs",
        "VLMC": "VLMCs",
        "ABC": "ABCs",
        "SMC": "SMCs",
        "OPC": "OPCs",
        "OEC": "OECs",
        "Ependymal":"Ependymal",
        "Astroependymal": "Astroependymal",
        "CHOR": "Choroid-Plexus",
        "Tanycyte": "Tanycytes",
        "Hypendymal": "Ependymal",
        "Monocytes": "Immune-Other",
        "DC": "Immune-Other",
        "Lymphoid": "Immune-Other",
        "BAM": "BAMs",
        "IMN": "Neurons-Granule-Immature",
        "Neurons-Glut-Chol": "Neurons-Other",
        "Neurons-Gaba-Glut": "Neurons-Other",
        "Neurons-Chol": "Neurons-Other",
        "Neurons-Hist-Gaba": "Neurons-Other",
        "Neurons-Gaba-Chol": "Neurons-Other",
        "Neurons-Glut-Sero": "Neurons-Other",
        "Neurons-Dopa-Gaba": "Neurons-Dopa",
        "Neurons-Gly-Gaba": "Neurons-Glyc-Gaba",
        "Neurons-Gaba-Glut": "Neurons-Gaba",
        "Neurons-Glut-Sero": "Neurons-Glut",
    }

    result = result.replace(replacement_dict)

    return result


def map_gene_symbols_to_ensembl(adata, species='mouse', logger=None):
    """
    Maps gene symbols in an AnnData object to Ensembl IDs.
    
    Args:
        adata: AnnData object with gene symbols in var['gene']
        species: Species name for querying (default: 'mouse')
        logger: Logger object (default: None)
        
    Returns:
        Updated AnnData object with Ensembl IDs in var['ensmus_id']
    """
    # load on demand
    import mygene
    
    mg = mygene.MyGeneInfo()
    
    # Initial mapping of gene symbols to Ensembl IDs
    gene_symbols = adata.var["gene"].unique().tolist()
    results = mg.querymany(gene_symbols, scopes='symbol', species=species, 
                          fields='ensembl.gene', returnall=True)
    
    gene_to_ensembl = {}
    for hit in results['out']:
        if hit.get('ensembl') and 'gene' in hit['ensembl']:
            gene_to_ensembl[hit['query']] = hit['ensembl']['gene']
        else:
            gene_to_ensembl[hit['query']] = 'NA'
    
    adata.var['ensmus_id'] = adata.var['gene'].map(gene_to_ensembl)
    
    # Handle missing Ensembl IDs
    missing_genes = adata.var[adata.var['ensmus_id'] == 'NA']['gene'].unique().tolist()
    
    # Query MyGene for extended info on missing genes
    missing_results = mg.querymany(missing_genes, scopes=['symbol', 'alias'], 
                                  species=species, fields=['ensembl.gene', 'symbol'])
    
    ensembl_replacements = {}
    symbol_replacements = {}
    for result in missing_results:
        if 'ensembl' in result and 'gene' in result['ensembl']:
            ensembl_replacements[result['query']] = result['ensembl']['gene']
            symbol_replacements[result['query']] = result['symbol']
        
        # If still not found, try expansive search
        if result['query'] not in ensembl_replacements:
            expanded_search = mg.query(result['query'] + "*", species=species, 
                                    fields=['ensembl.gene', 'symbol', 'alias'])
            
            if expanded_search['hits']:
                for hit in expanded_search['hits']:
                    aliases = hit.get('alias', [])
                    if isinstance(aliases, str):
                        aliases = [aliases]
                    
                    if result['query'] in aliases and 'ensembl' in hit and 'gene' in hit['ensembl']:
                        ensembl_replacements[result['query']] = hit['ensembl']['gene']
                        symbol_replacements[result['query']] = hit['symbol']
                        break
    
    # Apply replacements to both columns
    for old_gene, new_gene in symbol_replacements.items():
        mask = adata.var['gene'] == old_gene
        adata.var.loc[mask, 'gene'] = new_gene
        adata.var.loc[mask, 'ensmus_id'] = ensembl_replacements[old_gene]
    
    # Summary statistics
    genes_updated = len(symbol_replacements)
    remaining_na = len(adata.var[adata.var['ensmus_id'] == 'NA']['gene'].unique())

    logger.info("Genes updated: {}".format(genes_updated))
    logger.info("Remaining gene symbols: {}".format(remaining_na))
    
    return adata


def mark_low_quality_mappings(metadata, target_column, mad_factor, level):
    """Add new column for low-quality mappings marked as "Undefined" based on correlation MAD threshold as suggested by Allen Institute."""
    level_key = f"{target_column}_{level}"
    cor_key = f"{target_column}_cor_{level}"
    out_key = f"{target_column}_{level}_incl_low_quality"

    # Calculate the low-quality mask
    low_quality_mask = metadata.groupby(level_key)[cor_key].transform(
        lambda x: x < (x.median() - mad_factor * median_abs_deviation(x, nan_policy='omit'))
    )

    # Fill NaN values with False
    low_quality_mask = low_quality_mask.fillna(False)

    # Create new column with either "Undefined" or original value
    metadata[out_key] = metadata[level_key]
    metadata.loc[low_quality_mask, out_key] = "Undefined"
    

def normalize_counts(adata, method="area", target_sum=250, logger=None):
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
            logger.warning("Area not found, using volume instead.")
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
    adata.layers[norm_layer_name] = sp.csr_matrix(
        sp.csr_matrix(adata.X, dtype=np.float32).multiply(1 / size_factor[:, None]),
        dtype=np.float32,
    )

    # Scale to target sum
    norm_matrix = adata.layers[norm_layer_name]
    scaling_factors = target_sum / norm_matrix.sum(axis=1).A1
    adata.layers[norm_layer_name] = norm_matrix.multiply(scaling_factors[:, None])

    # Verify mean count matches target
    assert round(adata.layers[norm_layer_name].sum(axis=1).mean()) == target_sum, (
        f"normalized to {round(adata.layers[norm_layer_name].sum(axis=1).mean())}, norm_layer dtype: {type(adata.layers[norm_layer_name])}, adata shape: {adata.layers[norm_layer_name].shape}"
    )

    # Log-transform and z-score
    log_layer_name = f"{norm_layer_name.split('_')[0]}_log1p_norm"
    adata.layers[norm_layer_name] = sp.csr_matrix(adata.layers[norm_layer_name])
    adata.layers[log_layer_name] = sc.pp.log1p(adata.layers[norm_layer_name], copy=True)
    adata.layers["zscore"] = sc.pp.scale(
        adata.layers[log_layer_name], zero_center=True, max_value=None, copy=True
    )

    return adata


def plot_mad_thresholds(
    Allen_MMC_metadata,
    out_path,
    name="mad_threshold",
    group_column="allen_CLUS",
    value_column="allen_cor_CLUS",
    mad_factor=3,
    figsize=(13, 7),
):
    """Plot a violin plot of data grouped by a specified column, with horizontal lines showing MAD thresholds.

    Parameters:
    - Allen_MMC_metadata: DataFrame containing the data
    - group_column: Column name for grouping (e.g. "allen_CLUS")
    - value_column: Column name for the values (e.g. "allen_avg_cor_CLUS")
    - mad_factor: The factor for excluding outliers based on MAD (default is 3)
    """
    # Compute MAD within each cluster
    grouped_mad = Allen_MMC_metadata.groupby(group_column)[value_column].apply(
        median_abs_deviation
    )
    unique_groups = Allen_MMC_metadata[group_column].unique()

    with rc_context({"figure.figsize": figsize}):
        plt.grid(True, zorder=0)

        # Create violin plot
        sns.violinplot(
            data=Allen_MMC_metadata,
            x=group_column,
            y=value_column,
            density_norm="width",
            inner="quart",
            zorder=2,
        )

        # Add horizontal lines for each group's MAD threshold
        for group, mad in grouped_mad.items():
            group_median = Allen_MMC_metadata[
                Allen_MMC_metadata[group_column] == group
            ][value_column].median()
            threshold = group_median - mad_factor * mad

            plt.axhline(
                y=threshold,
                color="r",
                linestyle="-",
                xmin=(list(unique_groups).index(group) - 0) / len(unique_groups),
                xmax=(list(unique_groups).index(group) + 0.9) / len(unique_groups),
                zorder=2,
            )

        plt.xticks(rotation=90, ha="center")
        plt.tight_layout()
        plt.savefig(os.path.join(out_path, name + ".png"))
        plt.close()

        return


def process_adata(adata, logger):
    """
    Preprocesses AnnData object:
    - Filter cells
    - Normalize count data by cell volume
    - Compute PCA, neighbors, and UMAP
    - Copy Allen cell type annotations from obsm to obs
    """

    # Filter low-quality cells
    sc.pp.filter_cells(adata, min_counts=20)
    sc.pp.filter_cells(adata, min_genes=5)

    # Check for presence of cell type mapping
    assert not adata.obsm["allen_cell_type_mapping"].isna().any().all(), "All entries in 'allen_cell_type_mapping' are NaN"

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()

    # Normalize counts
    if "area" in adata.obs.columns:
        normalize_counts(adata, method="area")
    elif "volume" in adata.obs.columns:
        normalize_counts(adata, method="volume")

    normalize_counts(adata, method="library")

    selected_layer = "zscore"
    adata.X = adata.layers[selected_layer]
    logger.info(f"Using layer '{selected_layer}' for downstream analysis.")

    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Copy cell type mapping to obs
    prefix = "cell_type_mmc"
    cellmapping = {
        f"{prefix}_raw": "allen_SUBC",
        f"{prefix}_incl_low_quality": "allen_SUBC_incl_low_quality",
        f"{prefix}_is_mixed": "allen_SUBC_is_mixed",
        f"{prefix}_incl_mixed": "allen_SUBC_incl_mixed",
        f"{prefix}_mixed_names": "allen_SUBC_mixed_names",
        f"{prefix}_runner_up_1": "allen_runner_up_1_SUBC",
        f"{prefix}_runner_up_2": "allen_runner_up_2_SUBC",
        f"{prefix}_runner_up_1_incl_low_quality": "allen_runner_up_1_SUBC_incl_low_quality",
        f"{prefix}_runner_up_2_incl_low_quality": "allen_runner_up_2_SUBC_incl_low_quality",
        f"{prefix}_rup1_diff_prob": "allen_diff_rup1_prob_SUBC",
        f"{prefix}_rup2_diff_prob": "allen_diff_rup2_prob_SUBC",
    }

    for obs_key, allen_key in cellmapping.items():
        adata.obs[obs_key] = adata.obsm["allen_cell_type_mapping"].loc[adata.obs.index, allen_key]

    return adata


def plot_metric_distributions(df, out_path, file_name):
   """Plots histograms of raw MapMyCell metrics (correlations, probabilities, and runner-up metrics) for "CLAS" and "SUBC".
   
   For detailed explanation of MapMyCells output see https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/output.md."""
   # Sample data
   sample_data = df.sample(min(50000, len(df)))
   
   # Get columns ending with "CLAS" or "SUBC"
   cols_to_plot = [col for col in sample_data.select_dtypes(include=[np.number]).columns 
                   if col.endswith(("CLAS", "SUBC"))]
   
   # Set up plot grid
   n_cols, n_rows = 4, (len(cols_to_plot) + 3) // 4
   fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 8))
   axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
   
   # Plot histograms
   for i, col in enumerate(cols_to_plot):
       if i < len(axes):
           axes[i].hist(sample_data[col].dropna(), bins=80, range=(0, 1))
           axes[i].set_title(col, fontsize=8)
           axes[i].tick_params(labelsize=7)
   
   # Hide unused subplots
   for i in range(len(cols_to_plot), len(axes)):
       axes[i].set_visible(False)
   
   plt.tight_layout()
   plt.savefig(os.path.join(out_path, file_name + ".png"))
   plt.close()


def process_mapmycells_output(json_results):
    """Process MapMyCells JSON output into pd.DataFrame per cell with cell type labels, correlations, probabilities, and runner-up metrics.
    
    For detailed explanation of MapMyCells output see https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/output.md."""
    node_labels = {
        level: {
            node: json_results["taxonomy_tree"]["name_mapper"][level][node]["name"]
            for node in json_results["taxonomy_tree"][level]
        }
        for level in json_results["taxonomy_tree"]["hierarchy"]
    }
    cell_data = {c["cell_id"]: c for c in json_results["results"]}
    results = {}
    
    for level in json_results["taxonomy_tree"]["hierarchy"]:
        level_name = level.replace("CCN20230722_", "")
        for cell_id, cell_info in cell_data.items():
            if cell_id not in results:
                results[cell_id] = {}
            if level in cell_info:
                assignment = cell_info[level]["assignment"]
                prob = cell_info[level].get("bootstrapping_probability", np.nan)
                rups = cell_info[level].get("runner_up_assignment", [])
                rups_cor = cell_info[level].get("runner_up_correlation", [])
                rups_prob = cell_info[level].get("runner_up_probability", [])

                results[cell_id].update({
                    f"allen_{level_name}": node_labels[level].get(assignment, "Unknown"),
                    f"allen_cor_{level_name}": cell_info[level]["avg_correlation"],
                    f"allen_prob_{level_name}": prob,
                    f"allen_runner_up_1_{level_name}": node_labels[level].get(rups[0], "Unknown") if rups else np.nan,
                    f"allen_runner_up_1_cor_{level_name}": rups_cor[0] if rups_cor else np.nan,
                    f"allen_runner_up_1_prob_{level_name}": rups_prob[0] if rups_prob else np.nan,
                    f"allen_runner_up_2_{level_name}": node_labels[level].get(rups[1], "Unknown") if len(rups) > 1 else np.nan,
                    f"allen_runner_up_2_cor_{level_name}": rups_cor[1] if len(rups_cor) > 1 else np.nan,
                    f"allen_runner_up_2_prob_{level_name}": rups_prob[1] if len(rups_prob) > 1 else np.nan,
                })

                # Add columns with difference between main and follow-up probabilities, setting to 1 if any follow-up probability is NaN
                for i in [1, 2]:
                    rup_prob = results[cell_id].get(f"allen_runner_up_{i}_prob_{level_name}", np.nan)
                    diff_col = f"allen_diff_rup{i}_prob_{level_name}"
                    results[cell_id][diff_col] = prob - rup_prob if pd.notna(prob) and pd.notna(rup_prob) else 1.0

    return pd.DataFrame.from_dict(results, orient="index")


def revise_annotations(
    adata, 
    leiden_res=10.0,
    leiden_col=None, 
    cell_type_colors=None, 
    mmc_to_score_dict=None,
    score_threshold=0.5,
    top_n_genes=50,
    ABCAtlas_marker_df_path=None,
    logger=None
):
    """
    De-noise MapMyCells annotations by assigning cell types to Leiden clusters based on majority vote, plus revise annotations based on marker gene expression. Wrapper for score_cell_types(), assign_cell_types_to_clusters(), assign_final_cell_types()
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with gene expression and metadata.
    leiden_col : string
        Column containing leiden clusters
    cell_type_colors : dict
        Dictionary mapping cell types to colors.
    mmc_to_score_dict : dict
        Dictionary mapping cell type identifiers to their score cell types.
    score_threshold : float, default=0.5
        Threshold for cell type score to be considered valid.
    top_n_genes : int, default=50
        Number of top marker genes to use for scoring.

    """
    # Check if leiden clustering exists, run if not present
    if leiden_col not in adata.obs.columns:
        logger.info(f"Running Leiden clustering with resolution {leiden_res}")
        sc.tl.leiden(adata, key_added=leiden_col, resolution=leiden_res)
    else:
        logger.info(f"Using available Leiden clustering with resolution {leiden_res}")

    logger.info("\n2. Scoring marker gene expression...")
    
    # load and format markers
    ABCAtlas_marker_df = pd.read_csv(ABCAtlas_marker_df_path)
    cell_types = ABCAtlas_marker_df.columns.tolist()
    cell_type_dict = {}
    for cell_type in cell_types:
            if cell_type == "0":
                continue
            genes = ABCAtlas_marker_df[cell_type].iloc[0:].tolist()
            genes = [gene for gene in genes if pd.notna(gene)]
            cell_type_dict[cell_type] = genes
    del cell_type_dict["Bergmann"] # too few cells
        
    # Score cell types using marker genes
    adata = score_cell_types(
        adata, 
        marker_genes_dict=cell_type_dict,
        top_n_genes=top_n_genes
    )
    
    # Process each cell type annotation key
    annotation_keys = ["cell_type_mmc_raw", "cell_type_mmc_incl_mixed", "cell_type_mmc_incl_low_quality"]
    annotation_results = {}
    
    for key in annotation_keys:
        logger.info(f"Processing {key}")

        logger.info("3. Assigning cell types to Leiden clusters using majority vote")
        assigned_cell_type_dict, normalized_percentage = assign_cell_types_to_clusters(
            adata, 
            leiden_col=leiden_col, 
            cell_type_col=key
        )
        
        # Add to adata.obs
        adata.obs[f"{key}_clusters"] = adata.obs[leiden_col].map(assigned_cell_type_dict)

        logger.info("4. Revise majority vote using marker genes scores and identify final cell type")
        adata, undefined_reasons, cluster_score_matrix = assign_final_cell_types(
            adata, 
            cluster_labels_dict=assigned_cell_type_dict, 
            match_dict=mmc_to_score_dict, 
            leiden_col=leiden_col, 
            out_col=f"{key}_revised",
            score_threshold=score_threshold
        )

        # Print undefined reasons
        logger.info("Clusters marked as Undefined:")
        for cluster, reason in undefined_reasons.items():
            logger.info(f"Cluster {cluster}: {reason}")

        # Calculate summary statistics
        defined_count = (adata.obs[f"{key}_revised"] != "Undefined").sum()
        total_count = len(adata.obs)
        defined_percent = defined_count / total_count * 100

        logger.info(f"Annotation summary: {defined_count}/{total_count} cells ({defined_percent:.1f}%) assigned to cell types")

        # Update categorical values if colors are provided
        if cell_type_colors is not None:
            adata.obs[f"{key}_revised"] = pd.Categorical(
                adata.obs[f"{key}_revised"], categories=list(cell_type_colors.keys())
            )
            adata.obs[f"{key}_revised"] = adata.obs[f"{key}_revised"].cat.remove_unused_categories()

        counts = adata.obs[f"{key}_revised"].value_counts()
        logger.info(f"Value counts: {counts}")
        
        # Store results
        annotation_results[key] = {
            'defined_count': defined_count,
            'total_count': total_count,
            'defined_percent': defined_percent,
            'undefined_reasons': undefined_reasons
        }
    
    return adata, annotation_results, normalized_percentage


def run_mapmycells(adata, sample_name, method_name, annotation_path, data_dir, today):
    today = date.today().strftime("%Y%m%d")

    # Prepare adata
    adata_temp = adata.copy()
    adata_temp.var["gene"] = adata_temp.var.index
    adata_temp.var.index = adata_temp.var["ensmus_id"]
    if "counts" in adata_temp.layers:
        adata_temp.X = adata_temp.layers["counts"]
    elif len(adata_temp.layers) > 0:
        raise ValueError("'counts' layer not found, and other layers are present.")

    filename = f"{today}_adata_temp_mapmycells.h5ad"
    temp_file_path = os.path.join(annotation_path, filename)
    adata_temp.write(temp_file_path)

    # Ensure output directory exists
    output_dir = os.path.join(annotation_path, "mapmycells_out")
    os.makedirs(output_dir, exist_ok=True)

    # Build command
    ref_path = os.path.join(data_dir, "misc", "scRNAseq_ref_ABCAtlas_Yao2023Nature", "ref_taxonomy", "mouse_brain")
    cmd = [
        sys.executable, "-m", "cell_type_mapper.cli.from_specified_markers",
        "--query_path", os.path.join(annotation_path, filename),
        "--extended_result_path", os.path.join(output_dir, f"{today}_MapMyCells_{sample_name}_{method_name}.json"),
        "--csv_result_path", os.path.join(output_dir, f"{today}_MapMyCells_{sample_name}_{method_name}.csv"),
        "--drop_level", "CCN20230722_SUPT",
        "--cloud_safe", "False",
        "--query_markers.serialized_lookup", os.path.join(ref_path, "mouse_markers_230821.json"),
        "--precomputed_stats.path", os.path.join(ref_path, "precomputed_stats_ABC_revision_230821.h5"),
        "--type_assignment.normalization", "raw",
        "--type_assignment.n_processors", "4"
    ]

    subprocess.run(cmd, check=True)
    os.remove(temp_file_path)


def score_cell_types(adata, marker_genes_dict, top_n_genes=25, logger=None):
    """Score cells for marker gene expression using top n genes from marker_genes_dict."""
    if not marker_genes_dict:
        logger.warning("Empty marker_genes_dict. No scoring performed.")
        return adata
        
    for cell_type, genes in marker_genes_dict.items():
        if not genes:
            logger.warning(f"Empty gene list for {cell_type}. Skipping.")
            continue
            
        genes_to_use = genes[:top_n_genes] if len(genes) > top_n_genes else genes
        available_genes = [gene for gene in genes_to_use if gene in adata.var_names]
        
        if available_genes:
            score_name = f"score_{cell_type}"
            logger.info(f"Scoring {cell_type} with {len(available_genes)} genes")
            sc.tl.score_genes(adata, available_genes, score_name=score_name)
        else:
            logger.warning(f"No marker genes for {cell_type} found in dataset. Skipping.")
    
    return adata


def update_explorer(path, adata):
    """Update explorer files.

    Args:
    path: path to sdata.zarr and sdata.explorer
    adata: annotated adata
    """    

    # load on demand
    from spatialdata import read_zarr
    from sopa.io.explorer import write_cell_categories
    
    adata_original = read_zarr(os.path.join(path, "sdata.zarr"))["table"]
    adata_original.obs = adata.obs.merge(
        adata_original.obs, left_index=True, right_index=True, how="right"
    )
    write_cell_categories(os.path.join(path, "sdata.explorer"), adata_original)


