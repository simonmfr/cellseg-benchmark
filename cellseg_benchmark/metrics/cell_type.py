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


def plot_cell_type_distribution(cohort, results_suffix, methods=None, show=False):
    """Plot cell type distribution as stacked barplot, methods ordered by composition similarity.

    Uses celltype distribution of all samples together; "merged" averages cohorts equally and marks methods missing in a cohort with *.
    methods: optional list of clean method names to keep.
    """
    from scipy.cluster.hierarchy import leaves_list, linkage

    results_file = (
        Path(_constants.BASE_PATH)
        / "metrics"
        / cohort
        / "cell_type_metrics"
        / f"cell_type_distribution_{results_suffix}.csv"
    )
    cohorts = ["aging", "ABCAtlas", "VizgenMouseBrain"] if cohort == "merged" else [cohort]
    results_df = pd.concat(
        [pd.read_csv(str(results_file).replace(f"/{cohort}/", f"/{c}/"), index_col=0).assign(cohort=c) for c in cohorts], ignore_index=True
    ).query("sample == 'all'").drop(columns="sample")
    results_df["method"] = results_df.method.map(clean_method_name)
    if methods is not None:
        results_df = results_df[results_df.method.isin(methods)]
    partial = results_df.groupby("method").cohort.nunique() < len(cohorts)
    results_df = results_df.drop(columns="cohort").fillna(0).groupby("method").mean()
    results_df.index = [f"{m}*" if partial[m] else m for m in results_df.index]

    df_pct = results_df.T * 100
    df_pct = df_pct.reindex([ct for ct in _constants.cell_type_colors if ct in df_pct.index])
    rare = (df_pct.max(axis=1) < 0.5) & (df_pct.index != "Undefined")
    df_pct.loc["Other"] = df_pct[rare].sum()
    df_pct = df_pct[~rare.reindex(df_pct.index, fill_value=False)]
    colors = _constants.cell_type_colors | {"Other": "#9e9e9e"}
    df_pct = df_pct.iloc[:, leaves_list(linkage(df_pct.T, "average", optimal_ordering=True))][::-1]

    with plt.rc_context({"font.size": 7, "axes.linewidth": 0.6, "ytick.major.width": 0.6}):
        fig, ax = plt.subplots(figsize=(0.17 * df_pct.shape[1], 2.4))
        df_pct.T.plot.bar(stacked=True, ax=ax, color=[colors[ct] for ct in df_pct.index], width=0.82, lw=0, legend=False)
        ax.set(xlabel="", ylabel="Cells (%)", ylim=(0, 100), yticks=[0, 50, 100], xlim=(-0.6, df_pct.shape[1] - 0.4))
        ax.tick_params(axis="x", length=0, rotation=90)
        ax.spines[["top", "right", "bottom"]].set_visible(False)
        h, l = ax.get_legend_handles_labels()
        ax.legend(h[::-1], l[::-1], loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
                  handlelength=0.8, handleheight=0.8, labelspacing=0.2, borderaxespad=0)

    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(plot_path / f"cell_type_distribution_{results_suffix}{'_main' if methods else ''}.{ext}", dpi=300, bbox_inches="tight", transparent=True)
    if show:
        plt.show()
    plt.close(fig)
