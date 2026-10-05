import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd

from .. import _constants
from .utils import clean_method_name


def compute_cell_type_distribution(adata, celltype_name, **kwargs):
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


def compute_annotation_qc(adata, method, base_path=_constants.BASE_PATH, mad_factor=3.0, **kwargs):
    """`annotation_qc_summary` on post-merge cells, per sample and pooled ("all").

    Rebuilds per-cell MapMyCells output and SCALPEL QC from each sample's annotation folder,
    as in cell_type_annotation.py, and restricts it to the cells kept by merge_adata.
    """
    from ..cell_annotation_utils import (
        annotation_qc_summary,
        group_cell_types,
        process_mapmycells_output,
        scalpel_qc,
    )

    subs = {}
    for sample in adata.obs["sample"].unique():
        ann_path = Path(base_path) / "samples" / sample / "results" / method / "cell_type_annotation"
        jsons = sorted(ann_path.glob(f"mapmycells_out/*MapMyCells_{sample}_{method}.json"))
        if not jsons:
            print(f"No MapMyCells output for {sample}, skipping")
            continue
        with open(jsons[-1], "rb") as f:
            mmc = process_mapmycells_output(json.load(f))
        mmc = mmc.join(scalpel_qc(mmc["allen_cor_SUPT"], mmc["allen_SUPT"], mad_factor))
        mmc["allen_SUBC"] = group_cell_types(mmc["allen_SUBC"]).fillna("Undefined")
        mmc.index = mmc.index.astype(str)
        labels = pd.read_csv(
            ann_path / "adata_obs_annotated.csv",
            usecols=["cell_id", "cell_type_vote", "cell_type_revised"],
            dtype={"cell_id": str},
        ).set_index("cell_id")

        sub = adata[adata.obs["sample"] == sample]
        ids = sub.obs_names.astype(str)
        ids = ids.where(ids.isin(labels.index), ids.str.replace(r"-\d+$", "", regex=True))
        keep = ids.isin(labels.index) & ids.isin(mmc.index)
        if not keep.all():
            print(f"{sample}: {(~keep).sum()} merged cells without annotation, dropped")
        sub, ids = sub[keep], ids[keep]
        subs[sample] = ad.AnnData(
            obs=labels.loc[ids].set_axis(sub.obs_names).assign(volume=sub.obs["volume_final"].to_numpy()),
            obsm={"allen_cell_type_mapping": mmc.loc[ids].set_axis(sub.obs_names)},
            var=sub.var[[]],
            layers={"counts": sub.layers["counts"]},
        )
    if not subs:
        return None
    subs["all"] = ad.concat(subs.values(), merge="first")
    return pd.DataFrame({s: annotation_qc_summary(a) for s, a in subs.items()}).T.reset_index(names="sample")


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
        Path(_constants.BASE_PATH)
        / "metrics"
        / cohort
        / "cell_type_metrics"
        / f"cell_type_distribution_{results_suffix}.csv"
    )
    results_df = pd.read_csv(results_file, index_col=0)
    # clean method names for plotting
    results_df["method"] = results_df["method"].map(clean_method_name)
    results_df = (
        results_df[results_df["sample"] == "all"]
        .drop(columns=["sample"])
        .set_index("method")
    )

    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    df_pct = results_df.T * 100
    legend_order = list(_constants.cell_type_colors.keys())
    df_pct = df_pct.reindex([ct for ct in legend_order if ct in df_pct.index])
    if "Undefined" in df_pct.index:
        cols_sorted = df_pct.loc["Undefined"].sort_values(ascending=False).index
        df_pct = df_pct[cols_sorted]

    df_pct = df_pct[::-1]
    colors = [_constants.cell_type_colors[ct] for ct in df_pct.index]

    fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
    df_pct.T.plot(kind="bar", stacked=True, ax=ax, color=colors, width=0.8)

    ax.set_ylabel("% of Cells", fontsize=14)
    ax.set_xlabel("Segmentation Method", fontsize=14)
    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.yticks(fontsize=10)

    from matplotlib.patches import Patch

    legend_labels = df_pct[::-1].index.tolist()
    legend_colors = [_constants.cell_type_colors[ct] for ct in legend_labels]
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
