import ast
import warnings

import pandas as pd
import scanpy as sc
from shapely.geometry import Polygon


def compute_outlier_percentage(sdata, min_counts=25, min_genes=5):
    """Compute the percentage of outlier cells based on low counts/genes."""
    # min_counts=50 used in squidpy tutorial, 25 used in colab script from vizgen, 20 counts or 5 genes used in Allen 2023 paper
    outliers_percentage_dict = {}

    for name, table in sdata.tables.items():
        if name.startswith("adata_"):
            sc.pp.calculate_qc_metrics(
                table, qc_vars=[], percent_top=None, inplace=True
            )
            table.obs["cell_outlier"] = (table.obs["total_counts"] < min_counts) | (
                table.obs["n_genes_by_counts"] < min_genes
            )
            table.obs["cell_outlier"] = table.obs["cell_outlier"].astype(bool)
            outliers_count = table.obs["cell_outlier"].sum()
            outliers_percentage = (outliers_count / table.shape[0]) * 100
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
        name.replace("adata_", ""): table[~table.obs["cell_outlier"]]
        .obs["n_genes_by_counts"]
        .mean()
        for name, table in sdata.tables.items()
        if name.startswith("adata_")
    }


def mean_counts_per_cell(sdata):
    """Calculate the mean counts per cell for each table and return as dict.

    Requires output from compute_outlier_percentage().
    """
    return {
        name.replace("adata_", ""): table[~table.obs["cell_outlier"]]
        .obs["total_counts"]
        .mean()
        for name, table in sdata.tables.items()
        if name.startswith("adata_")
    }


def mean_cell_area_2d(sdata, zindex=3):
    """Compute mean 2D cell area per dataset, excluding outlier cells, using the specified ZIndex if available."""
    mean_area_dict = {}

    for name, boundaries in sdata.shapes.items():
        if name.startswith("boundaries_"):
            table_name = name.replace("boundaries_", "adata_")  # Match table name
            if table_name not in sdata.tables:
                warnings.warn(
                    f"Missing corresponding table {table_name} for {name}. Skipping dataset."
                )
                continue

            table = sdata.tables[table_name]
            if "cell_outlier" not in table.obs:
                warnings.warn(
                    f"Missing 'cell_outlier' column in {table_name}. Using all cells."
                )
                valid_cells = table.obs.index  # Use all cells if no outlier info
            else:
                valid_cells = table.obs.loc[
                    ~table.obs["cell_outlier"]
                ].index  # Exclude outliers

            if boundaries.empty or "geometry" not in boundaries:
                continue  # Skip if data is missing

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
                    warnings.warn(
                        f"Only one z-plane detected. Setting zindex to {zindex}. Verify sdata loading if incorrect."
                    )

                boundaries = boundaries[boundaries["ZIndex"] == zindex]
            else:
                warnings.warn(
                    f"ZIndex column missing in {name}. Using the entire dataset."
                )

            # Filter boundaries using index
            boundaries = boundaries.loc[boundaries.index.intersection(valid_cells)]

            # Compute mean cell area
            mean_area_dict[name.replace("boundaries_", "")] = (
                boundaries["geometry"].apply(lambda x: x.area).mean()
            )

    return mean_area_dict


def combine_metrics(sdata):
    """Run selected metrics and merge results into a dataframe."""
    # Get the individual metrics as dicts
    outliers_dict = compute_outlier_percentage(sdata)
    cell_counts_dict = count_cells(sdata)
    mean_genes_dict = mean_genes_per_cell(sdata)
    mean_counts_dict = mean_counts_per_cell(sdata)
    mean_cell_area_2d_dict = mean_cell_area_2d(sdata)

    # Convert each dict into a df and join them
    outliers_df = pd.DataFrame(
        list(outliers_dict.items()), columns=["Dataset", "Outlier Percentage"]
    ).set_index("Dataset")
    cell_counts_df = pd.DataFrame(
        list(cell_counts_dict.items()), columns=["Dataset", "Cell Count"]
    ).set_index("Dataset")
    mean_genes_df = pd.DataFrame(
        list(mean_genes_dict.items()), columns=["Dataset", "Mean Genes per Cell"]
    ).set_index("Dataset")
    mean_counts_df = pd.DataFrame(
        list(mean_counts_dict.items()), columns=["Dataset", "Mean Counts per Cell"]
    ).set_index("Dataset")
    mean_cell_area_2d_df = pd.DataFrame(
        list(mean_cell_area_2d_dict.items()), columns=["Dataset", "Mean Cell Area 2D"]
    ).set_index("Dataset")

    # Merge all dfs on the 'Dataset' index
    combined_df = outliers_df.join(
        [cell_counts_df, mean_genes_df, mean_counts_df, mean_cell_area_2d_df]
    )

    return combined_df
