import pandas as pd


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
