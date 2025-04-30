import ast
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from shapely.geometry import Polygon


def compute_outlier_percentage(sdata, min_counts=25, min_genes=5, inplace=True):
    """Compute percentage of outlier cells based on low counts/genes.

    If inplace = True, adds qc metrics and boolean "cell_outlier" to adata.obs. No filtering is done.

    Parameters: sdata (SpatialData), min_counts (int, default=25),
    min_genes (int, default=5), inplace (bool, default=True)
    Returns: dict of dataset names to outlier percentages
    """
    outliers_percentage_dict = {}

    for name, table in sdata.tables.items():
        if name.startswith("adata_"):
            if not inplace:
                table = table.copy()

            sc.pp.calculate_qc_metrics(
                table, qc_vars=[], percent_top=None, inplace=True
            )
            table.obs["cell_outlier"] = (table.obs["total_counts"] < min_counts) | (
                table.obs["n_genes_by_counts"] < min_genes
            )
            table.obs["cell_outlier"] = table.obs["cell_outlier"].astype(bool)
            outliers_count = table.obs["cell_outlier"].sum()
            outliers_percentage = round((outliers_count / table.shape[0]) * 100, 1)
            outliers_percentage_dict[name.replace("adata_", "")] = outliers_percentage

            if not inplace:
                sdata.tables[name] = table

    return outliers_percentage_dict


def count_cells(sdata):
    """Count non-outlier cells for each table in sdata and return as dict.

    Requires output from compute_outlier_percentage().
    """
    try:
        return {
            name.replace("adata_", ""): table[~table.obs["cell_outlier"]].shape[0]
            for name, table in sdata.tables.items()
            if name.startswith("adata_")
        }
    except Exception as e:
        warnings.warn(f"Error in count_cells: {str(e)}")
        return {}


def mean_genes_per_cell(sdata):
    """Calculate the mean genes per cell for each table and return as dict.

    Requires output from compute_outlier_percentage().
    """
    try:
        return {
            name.replace("adata_", ""): int(
                table[~table.obs["cell_outlier"]].obs["n_genes_by_counts"].mean()
            )
            for name, table in sdata.tables.items()
            if name.startswith("adata_")
        }
    except Exception as e:
        warnings.warn(f"Error in mean_genes_per_cell: {str(e)}")
        return {}


def mean_counts_per_cell(sdata):
    """Calculate the mean counts per cell for each table and return as dict.

    Requires output from compute_outlier_percentage().
    """
    try:
        return {
            name.replace("adata_", ""): int(
                table[~table.obs["cell_outlier"]].obs["total_counts"].mean()
            )
            for name, table in sdata.tables.items()
            if name.startswith("adata_")
        }
    except Exception as e:
        warnings.warn(f"Error in mean_counts_per_cell: {str(e)}")
        return {}


def compute_cell_areas_2d(sdata, zindex=3, column_name="area", verbose=True):
    """Compute and add 2D cell area to adata.obs for each dataset in SpatialData.

    Parameters:
    -----------
    sdata : SpatialData
        The SpatialData object containing shape and table data
    zindex : int, default=3
        Z-plane index to use if ZIndex column is present
    column_name : str, default="area"
        Name of the new column in adata.obs

    Returns:
    --------
    None
        Modifies adata.obs by adding the cell area column for each dataset
    """
    for boundaries_name, boundaries in sdata.shapes.items():
        if not boundaries_name.startswith("boundaries_"):
            continue

        dataset_name = boundaries_name.replace("boundaries_", "")
        table_name = f"adata_{dataset_name}"

        if table_name not in sdata.tables:
            warnings.warn(
                f"Missing corresponding adata {table_name} for {boundaries_name}. Skipping dataset."
            )
            continue

        adata = sdata.tables[table_name]

        if boundaries.empty or "geometry" not in boundaries:
            warnings.warn(
                f"Empty or invalid geometry data in {boundaries_name}. Skipping dataset."
            )
            continue

        # Convert geometry from string to Polygon if needed
        boundaries = boundaries.copy()
        boundaries["geometry"] = boundaries["geometry"].apply(
            lambda x: Polygon(ast.literal_eval(x)) if isinstance(x, str) else x
        )

        # Handle ZIndex if present
        if "ZIndex" in boundaries:
            unique_z = boundaries["ZIndex"].unique()
            if len(unique_z) == 1:
                zindex_0 = unique_z[0]
                print(
                    f"Only one z-plane detected for {dataset_name}. Using zindex {zindex_0}."
                )
            boundaries = boundaries[boundaries["ZIndex"] == zindex_0]

        # Compute area for each cell
        areas = boundaries["geometry"].apply(lambda x: x.area)

        # Add area to adata.obs, matching cells by index
        areas.index = areas.index.astype(str)
        common_indices = adata.obs.index.intersection(areas.index)
        adata.obs[column_name] = np.nan
        adata.obs.loc[common_indices, column_name] = areas.loc[common_indices]

        if verbose:
            if "ZIndex" in boundaries and len(unique_z) == 1:
                print(f"Added {column_name} for {dataset_name} (z={zindex_0})")
            else:
                print(f"Added {column_name} for {dataset_name} (z={zindex})")


def mean_cell_area_2d(sdata, column_name="area"):
    """Compute mean 2D cell area per dataset across cells, excluding outlier cells.

    Requires cell area in adata.obs, such as from running compute_cell_areas_2d.
    """
    mean_area_dict = {}

    for name, boundaries in sdata.shapes.items():
        if name.startswith("boundaries_"):
            dataset_name = name.replace("boundaries_", "")
            table_name = f"adata_{dataset_name}"
            if table_name not in sdata.tables:
                warnings.warn(f"Missing table {table_name} for {name}. Skipping.")
                continue

            adata = sdata.tables[table_name]
            if column_name not in adata.obs:
                warnings.warn(
                    f"{column_name} missing in {table_name}. Skipping. Did you run compute_cell_areas_2d first?"
                )
                continue

            # Filter valid cells (exclude outliers if column exists)
            if "cell_outlier" in adata.obs:
                valid_cells = adata.obs.loc[~adata.obs["cell_outlier"], column_name]
            else:
                warnings.warn(
                    f"'cell_outlier' column missing in {table_name}. Using all cells."
                )
                valid_cells = adata.obs[column_name]

            # Compute mean area, ignoring NaNs, and round to integer
            mean_value = valid_cells.mean()
            if pd.notna(mean_value):
                mean_area_dict[dataset_name] = round(mean_value)

    return mean_area_dict


def transcript_density(cell_volume_3d, n_transcripts):
    """Calculate the transcript density in 3D.

    Calculation done by dividing the number of transcripts
    per cell by the cell volume.

    Parameters:
        cell_volume_3d (pd.DataFrame): Containing the 3D volume of each cell.
        n_transcripts (pd.DataFrame): Containing the total number of transcripts per cell.

    Returns:
        pandas.Series: The transcript density (transcripts per unit volume) for each cell.
    """
    cell_volume_3d.index = cell_volume_3d.index.astype(str)
    n_transcripts.index = n_transcripts.index.astype(str)

    return n_transcripts.loc[cell_volume_3d.index] / cell_volume_3d


def mean_transcript_densities(sdata, area_col="area"):
    """Calculate transcript densities (transcripts/cell_volume_3d) for all tables in sdata.

    Requires area_col (default: "area") and "total_counts" in adata.obs.
    """
    transcript_densities = {}
    for table_name in sdata.tables.keys():
        # Strip "adata_" prefix from table name for consistency
        clean_name = table_name.replace("adata_", "")

        # Check required columns
        if area_col not in sdata[table_name].obs:
            warnings.warn(f"Missing {area_col} in {table_name}. Skipping.")
            continue
        if "total_counts" not in sdata[table_name].obs:
            warnings.warn(f"Missing 'total_counts' in {table_name}. Skipping.")
            continue

        # Filter valid cells (exclude outliers if column exists)
        if "cell_outlier" in sdata[table_name].obs:
            valid_cells = ~sdata[table_name].obs["cell_outlier"]
        else:
            warnings.warn(
                f"'cell_outlier' column missing in {table_name}. Using all cells."
            )
            valid_cells = pd.Series(True, index=sdata[table_name].obs.index)

        density = transcript_density(
            cell_volume_3d=sdata[table_name].obs.loc[valid_cells][area_col] * 7,
            n_transcripts=sdata[table_name].obs.loc[valid_cells]["total_counts"],
        )
        transcript_densities[clean_name] = round(density.mean(), 4)

    return transcript_densities


def combine_metrics(sdata):
    """Run and combine selected metrics into a single dataframe."""
    # Compute cell areas
    compute_cell_areas_2d(sdata, zindex=3, column_name="area", verbose=True)

    # Define metrics
    metric_dicts = {
        "% Low Quality Cells": compute_outlier_percentage(sdata),
        "Number of Cells": count_cells(sdata),
        "Genes per Cell (mean)": mean_genes_per_cell(sdata),
        "Counts per Cell (mean)": mean_counts_per_cell(sdata),
        "Cell Area (mean in z=3)": mean_cell_area_2d(
            sdata
        ),  # requires compute_cell_areas_2d()
        "Cell Transcript Density (mean per um3)": mean_transcript_densities(sdata),
    }

    # Get all unique dataset names
    all_datasets = sorted(
        set().union(*[dict.keys() for dict_i in metric_dicts.values()])
    )

    # Create DataFrame with datasets as index
    df = pd.DataFrame(index=all_datasets)
    df.index.name = "Dataset"

    # Fill each column with corresponding metric values
    for col_name, metric_dict in metric_dicts.items():
        df[col_name] = [metric_dict.get(dataset, pd.NA) for dataset in all_datasets]

    # Convert integer-like columns to Int64 type (nullable integer)
    for col in df.columns:
        try:
            if all(float(x).is_integer() for x in df[col] if pd.notna(x)):
                df[col] = df[col].astype("Int64")
        except (ValueError, AttributeError):
            pass

    return df
