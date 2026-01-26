from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from cellseg_benchmark._constants import BASE_PATH, cell_type_colors


def compute_cell_type_distribution(adata, celltype_name):
    """Compute distribution of celltypes per sample.

    Args:
        adata: anndata to compute score with
        celltype_name: name of celltype column in adata.obs.

    Returns:
        results DataFrame or None
    """
    # check if celltype name exists
    if celltype_name not in adata.obs.columns and (celltype_name != "leiden"):
        print(f"{celltype_name} not found in adata.obs")
        return None
    # compute clustering score per sample
    results = _cell_type_distribution(adata, celltype_name)
    return results.reset_index(names="sample")


def _cell_type_distribution(adata, celltype_name):
    results = pd.DataFrame(columns=adata.obs[celltype_name].unique())
    for sample in adata.obs["sample"].unique():
        cur_adata = adata[adata.obs["sample"] == sample]
        results.loc[sample] = dict(
            cur_adata.obs[celltype_name].value_counts()
            / len(cur_adata.obs[celltype_name])
        )
    # compute for all samples together
    results.loc["all"] = dict(
        adata.obs[celltype_name].value_counts() / len(adata.obs[celltype_name])
    )
    results = results.fillna(0)
    return results


def plot_cell_type_distribution(cohort, results_suffix, show=False):
    """Plot cell type distribution as stacked barplot.

    Uses celltype distribution of all samples together.
    """
    results_file = (
        Path(BASE_PATH)
        / "metrics"
        / cohort
        / "cell_type_metrics"
        / f"cell_type_distribution_{results_suffix}.csv"
    )
    results_df = pd.read_csv(results_file, index_col=0)
    results_df = (
        results_df[results_df["sample"] == "all"]
        .drop(columns=["sample"])
        .set_index("method")
    )

    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    df_pct = results_df.T * 100
    legend_order = list(cell_type_colors.keys())
    df_pct = df_pct.reindex([ct for ct in legend_order if ct in df_pct.index])
    if "Undefined" in df_pct.index:
        cols_sorted = df_pct.loc["Undefined"].sort_values(ascending=False).index
        df_pct = df_pct[cols_sorted]

    df_pct = df_pct[::-1]
    colors = [cell_type_colors[ct] for ct in df_pct.index]

    fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
    df_pct.T.plot(kind="bar", stacked=True, ax=ax, color=colors, width=0.8)

    ax.set_ylabel("% of Cells", fontsize=14)
    ax.set_xlabel("Segmentation Method", fontsize=14)
    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.yticks(fontsize=10)

    from matplotlib.patches import Patch

    legend_labels = df_pct[::-1].index.tolist()
    legend_colors = [cell_type_colors[ct] for ct in legend_labels]
    handles = [
        Patch(facecolor=color, label=label)
        for color, label in zip(legend_colors, legend_labels)
    ]

    ax.legend(
        handles=handles,
        title="Cell Type",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        fontsize=9,
    )
    plt.tight_layout()

    if show:
        plt.show()
    fig.savefig(plot_path / f"cell_type_distribution_{results_suffix}.png")
    plt.close(fig)
