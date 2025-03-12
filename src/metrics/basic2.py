import pandas as pd


def transcript_density_3d(cell_volume_3d, n_transcripts):
    """Calculate the transcript density in 3D.

    Calculation done by dividing the number of transcripts
    per cell by the cell volume.

    Parameters:
        cell_volume_3d (pd.DataFrame): C containing the 3D volume of each cell.
        n_transcripts (pd.DataFrame): Containing the total number of transcripts per cell.

    Returns:
        pandas.Series: The transcript density (transcripts per unit volume) for each cell.
    """
    cell_volume_3d.index = cell_volume_3d.index.astype(int)
    n_transcripts.index = n_transcripts.index.astype(int)

    return n_transcripts.loc[cell_volume_3d.index] / cell_volume_3d


def pct_unassigned_transcripts(
    transcripts_df, no_cell_key=-1, cell_identifier="cell_id"
):
    """Calculate the percentage of all transcripts not assigned to any cell.

    Parameters:
        transcripts_df (pd.DataFrame): DataFrame containing transcript data including cell assignment
        no_cell_key (int): The value representing no cell type match (default is -1).
        cell_identifier (str): The column name for cell assignment (default is "cell_id").

    Returns:
        float: Percentage of transcripts not assigned to the excluded cell type.
    """
    pct = (
        (transcripts_df[cell_identifier] == no_cell_key).sum() / len(transcripts_df)
    ) * 100
    return round(pct, 2)


def top_unsegmented_genes(
    transcripts_df, no_cell_key=-1, cell_identifier="cell_id", gene_identifier="gene"
):
    """Ranks genes based on the proportion of unsegmented vs segmented transcripts.

    Parameters:
        transcripts_df (pd.DataFrame): A dataframe containing transcript data with cell assignment and gene labels.
        no_cell_key (int, optional): A value indicating the absence of a segmented cell (default is -1).
        cell_identifier (str, optional): The column name for cell identifiers (default is "cell_id").
        gene_identifier (str, optional): The column name for gene identifiers (default is "gene").

    Returns:
        pandas.DataFrame: A dataframe with the counts of segmented and unsegmented occurrences for each gene,
            sorted by the proportion of unsegmented vs segmented descending order.
    """
    segmented_counts = transcripts_df[transcripts_df[cell_identifier] != no_cell_key][
        gene_identifier
    ].value_counts()
    unsegmented_counts = transcripts_df[transcripts_df[cell_identifier] == no_cell_key][
        gene_identifier
    ].value_counts()

    gene_counts = pd.DataFrame(
        {"segmented_count": segmented_counts, "unsegmented_count": unsegmented_counts}
    ).fillna(0)

    gene_counts["total_count"] = gene_counts.sum(axis=1)
    gene_counts["unsegmented_proportion"] = (
        gene_counts["unsegmented_count"] / gene_counts["total_count"]
    ).round(2)
    gene_counts.index.name = None

    return gene_counts.sort_values("unsegmented_proportion", ascending=False)
