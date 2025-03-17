import ast
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from shapely.geometry import Polygon


def compute_outlier_percentage(sdata, min_counts=25, min_genes=5):
    """Compute the percentage of outlier cells based on low counts/genes.

    Flags outliers in adata.obs["cell_outlier"], but does not filter.

    min_counts=50 used in squidpy tutorial, 25 used in colab script from vizgen, 20 counts or 5 genes used in Allen 2024.
    """
    outliers_percentage_dict = {}

    for name, table in sdata.tables.items():
        if name.startswith("adata_"):
            sc.pp.calculate_qc_metrics(
                table, qc_vars=[], percent_top=None, inplace=True
            )
            table.obs["cell_outlier"] = (table.obs["total_counts"] < min_counts) | (
                table.obs["n_genes_by_counts"] < min_genes
            )

            # additionally filter by volume? Done e.g. by Allen 2023 Cell
            # (if cell volume is < 100 um3 or >3 times the median across all cells (across all experiments))

            table.obs["cell_outlier"] = table.obs["cell_outlier"].astype(bool)
            outliers_count = table.obs["cell_outlier"].sum()
            outliers_percentage = round((outliers_count / table.shape[0]) * 100, 1)
            outliers_percentage_dict[name.replace("adata_", "")] = outliers_percentage

    return outliers_percentage_dict


def count_cells(sdata):
    """Count non-outlier cells  for each table and return as dict.

    Requires output from compute_outlier_percentage().
    """
    return {
        name.replace("adata_", ""): table[~table.obs["cell_outlier"]].shape[0]
        for name, table in sdata.tables.items()
        if name.startswith("adata_")
    }


def mean_genes_per_cell(sdata):
    """Calculate the mean genes per cell for each table and return as dict.

    Requires output from compute_outlier_percentage().
    """
    return {
        name.replace("adata_", ""): int(
            table[~table.obs["cell_outlier"]].obs["n_genes_by_counts"].mean()
        )
        for name, table in sdata.tables.items()
        if name.startswith("adata_")
    }


def mean_counts_per_cell(sdata):
    """Calculate the mean counts per cell for each table and return as dict.

    Requires output from compute_outlier_percentage().
    """
    return {
        name.replace("adata_", ""): int(
            table[~table.obs["cell_outlier"]].obs["total_counts"].mean()
        )
        for name, table in sdata.tables.items()
        if name.startswith("adata_")
    }


def compute_cell_areas_2d(sdata, zindex=3, column_name="cell_area", verbose=True):
    """Compute 2D cell area for all adatas in SpatialData and add it to respective adata.obs.

    Parameters:
    -----------
    sdata : SpatialData
        The SpatialData object containing shapes and tables
    zindex : int, default=3
        The z-plane to use if ZIndex column is available
    column_name : str, default="cell_area"
        Name of the column to add to adata.obs

    Returns:
    --------
    None
        Modifies adata.obs in-place by adding the cell area column for each dataset
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
                zindex = unique_z[0]
                print(
                    f"Only one z-plane detected for {dataset_name}. Using zindex {zindex}."
                )
            boundaries = boundaries[boundaries["ZIndex"] == zindex]

        # Compute area for each cell
        areas = boundaries["geometry"].apply(lambda x: x.area)

        # Add area to adata.obs, matching cells by index
        common_indices = adata.obs.index.intersection(areas.index)
        adata.obs[column_name] = np.nan
        adata.obs.loc[common_indices, column_name] = areas.loc[common_indices]

        if verbose:
            print(f"Added {column_name} for {dataset_name}")


def mean_cell_area_2d(sdata, column_name="cell_area_2d"):
    """Compute mean 2D cell area per dataset across cells, excluding outlier cells.

    Requires cell area in adata.obs, such as from running compute_cell_areas_2d.
    """
    mean_area_dict = {}

    for name, boundaries in sdata.shapes.items():
        if name.startswith("boundaries_"):
            dataset_name = name.replace("boundaries_", "")
            table_name = f"adata_{dataset_name}"

            if table_name not in sdata.tables:
                warnings.warn(
                    f"Missing corresponding table {table_name} for {name}. Skipping dataset."
                )
                continue

            adata = sdata.tables[table_name]

            if column_name not in adata.obs:
                warnings.warn(
                    f"{column_name} column missing in {table_name}. Skipping dataset."
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
            mean_area = round(valid_cells.mean())
            mean_area_dict[dataset_name] = mean_area

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


def mean_transcript_densities(sdata):
    """Calculate transcript densities for all tables in sdata.

    Requires "cell_area_2d" and "total_counts" in adata.obs.
    """
    transcript_densities = {}

    for table_name in sdata.tables.keys():
        # Filter valid cells (exclude outliers if column exists)
        if "cell_outlier" in sdata[table_name].obs:
            valid_cells = ~sdata[table_name].obs["cell_outlier"]
        else:
            warnings.warn(
                f"'cell_outlier' column missing in {table_name}. Using all cells."
            )
            valid_cells = pd.Series(True, index=sdata[table_name].obs.index)

        density = transcript_density(
            cell_volume_3d=sdata[table_name].obs.loc[valid_cells]["cell_area_2d"] * 7,
            n_transcripts=sdata[table_name].obs.loc[valid_cells]["total_counts"],
        )

        transcript_densities[table_name] = round(density.mean(), 4)

    return transcript_densities


def combine_metrics(sdata):
    """Run selected metrics and merge results into a dataframe."""
    # Get the individual metrics as dicts
    outliers_dict = compute_outlier_percentage(sdata)
    cell_counts_dict = count_cells(sdata)
    mean_genes_dict = mean_genes_per_cell(sdata)
    mean_counts_dict = mean_counts_per_cell(sdata)
    compute_cell_areas_2d(sdata, zindex=3, column_name="cell_area_2d", verbose=False)
    mean_cell_area_2d_dict = mean_cell_area_2d(sdata)
    transcript_densities_dict = mean_transcript_densities(sdata)

    # Convert each dict into a df and join them
    outliers_df = pd.DataFrame(
        list(outliers_dict.items()), columns=["Dataset", "% Low Quality Cells"]
    ).set_index("Dataset")
    cell_counts_df = pd.DataFrame(
        list(cell_counts_dict.items()), columns=["Dataset", "Number of Cells"]
    ).set_index("Dataset")
    mean_genes_df = pd.DataFrame(
        list(mean_genes_dict.items()), columns=["Dataset", "Genes per Cell (mean)"]
    ).set_index("Dataset")
    mean_counts_df = pd.DataFrame(
        list(mean_counts_dict.items()), columns=["Dataset", "Counts per Cell (mean)"]
    ).set_index("Dataset")
    mean_cell_area_2d_df = pd.DataFrame(
        list(mean_cell_area_2d_dict.items()), columns=["Dataset", "Cell Area 2D (mean)"]
    ).set_index("Dataset")
    mean_transcript_densities_df = pd.DataFrame(
        list(transcript_densities_dict.items()),
        columns=["Dataset", "Cell Transcript Density (mean per um3)"],
    ).set_index("Dataset")

    # Merge all dfs on the 'Dataset' index
    combined_df = outliers_df.join(
        [
            cell_counts_df,
            mean_genes_df,
            mean_counts_df,
            mean_cell_area_2d_df,
            mean_transcript_densities_df,
        ]
    )

    return combined_df
