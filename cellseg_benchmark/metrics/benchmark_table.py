"""Benchmark summary tables: score all methods per cohort and draw them as funky heatmaps."""

from pathlib import Path

import numpy as np
import pandas as pd

from .. import _constants
from .utils import clean_method_name

categories = {
    "clustering": ["silhouette_score"],
    "marker": ["MECR_adj", "marker_f1"],
    "ficture": ["ficture_f1_score"],
    "z": ["Ovrlpy_stats_mean_integrity"],
}
descriptive = [
    "volume_final",
    "elongation",
    "intensities_DAPI",
    "cells_per_mm2",
    "pct_assigned_qced",
]
metrics = sum(categories.values(), []) + descriptive
exclude = ["Negative_Control_Visium", "Ficture_segments"]
metric_tables = {
    "general_stats": "general_stats",
    "assigned_transcript_counts": "assigned_transcripts",
    "cell_density": "cell_density",
    "MECR_score_all": "marker_gene_metrics",
    "ficture_f1": "ficture",
    "marker_f1_score_cell_type_revised": "marker_gene_metrics",
    "clustering_score_adata_integrated_cell_type_revised": "cell_type_metrics",
}
group_tables = {
    "MECR_score_groups": "marker_gene_metrics",
    "clustering_score_groups": "cell_type_metrics",
}


def add_scores(df):
    """Add total and expression scores as means of category means, sorted by total."""
    category_scores = pd.DataFrame(
        {name: df[cols].mean(axis=1) for name, cols in categories.items()}
    )
    df["total_score"] = category_scores.mean(axis=1)
    df["expr_score"] = category_scores[["clustering", "marker"]].mean(axis=1)
    return df.sort_values("total_score", ascending=False)


def build_benchmark_table(dfs, group=None, ref=None):
    """Average metrics per method and min-max scale them across non-control methods.

    Args:
        dfs: metric tables as loaded by build_benchmark_tables.
        group: key of _constants.celltype_groups to restrict the metrics to, or None for all cells.
        ref: (MECR fit, min, max) to scale with instead of fitting on this table.

    Returns:
        The scored table and its (MECR fit, min, max).
    """
    general = dfs["general_stats"]
    stats = [
        "volume_final",
        "elongation",
        "intensities_DAPI",
        "Ovrlpy_stats_mean_integrity",
    ]
    all_cells = general.query("cell_type_revised == 'all'").set_index(
        ["sample", "method"]
    )
    marker = dfs["marker_f1_score_cell_type_revised"]
    ficture = dfs["ficture_f1"][~dfs["ficture_f1"].cell_type.str.startswith("f1_")]
    if group is None:
        cells = all_cells
        transcripts = None
        mecr = dfs["MECR_score_all"]
        silhouette = dfs["clustering_score_adata_integrated_cell_type_revised"]
    else:
        types = _constants.celltype_groups[group]
        rows = general[general.cell_type_revised.isin(types)]
        stats_w = stats + ["n_counts"]
        cells = (
            rows[stats_w]
            .mul(rows.n_cells, axis=0)
            .assign(n_cells=rows.n_cells, sample=rows["sample"], method=rows.method)
            .groupby(["sample", "method"])
            .sum(min_count=1)
        )
        cells[stats_w] = cells[stats_w].div(cells.n_cells, axis=0)
        transcripts = cells.n_counts * cells.n_cells
        mecr = dfs["MECR_score_groups"].query("group == @group")
        silhouette = dfs["clustering_score_groups"].query("group == @group")
        marker = marker[
            marker.cell_type.isin(
                {_constants.merged_celltypes.get(t, t) for t in types}
            )
        ]
        ficture = ficture[ficture.cell_type.isin(types)]

    counts = (
        dfs["assigned_transcript_counts"]
        .groupby(["sample", "method"])[["assigned_count_qced", "total_count"]]
        .sum()
    )
    if transcripts is None:
        transcripts = counts.assigned_count_qced
    counts["tpc"] = transcripts / cells.n_cells
    counts["pct_assigned_qced"] = transcripts / counts.total_count
    density = dfs["cell_density"].groupby("method")[["n_cells", "tissue_mm2"]].sum()
    fraction = (
        cells.n_cells.groupby(level="method").sum()
        / all_cells.n_cells.groupby(level="method").sum()
    )

    df = pd.concat(
        [
            cells.groupby(level="method")[stats].mean(),
            counts.groupby(level="method")[["pct_assigned_qced", "tpc"]].mean(),
            silhouette.groupby("method")[["silhouette_score"]].mean(),
            mecr.query("gene1 == 'all'").groupby("method")[["MECR"]].mean(),
            marker.groupby(["sample", "method", "cell_type"])
            .f1_score.mean()
            .groupby(["sample", "method"])
            .mean()
            .groupby("method")
            .mean()
            .rename("marker_f1"),
            ficture.groupby(["sample", "method"])
            .f1.mean()
            .groupby("method")
            .mean()
            .rename("ficture_f1_score"),
            (density.n_cells / density.tissue_mm2 * fraction).rename("cells_per_mm2"),
        ],
        axis=1,
    ).drop(index=exclude, errors="ignore")

    df["is_control"] = df.index.str.startswith("Negative_Control")
    if ref is None:
        fit_data = df[~df.is_control].dropna(subset=["MECR", "tpc"])
        fit = np.polyfit(np.log10(fit_data.tpc), fit_data.MECR, 1)
    else:
        fit, low, high = ref
    df["MECR_adj"] = df.MECR - np.polyval(fit, np.log10(df.tpc))

    is_3d = df.index.str.startswith(tuple(_constants.methods_3D))
    is_ficture_seg = df.index.str.startswith("Ficture_segments")
    df.loc[is_ficture_seg, "ficture_f1_score"] = np.nan
    df["na_z"] = is_3d & df.Ovrlpy_stats_mean_integrity.isna()
    df["na_ficture"] = (is_3d | is_ficture_seg) & df.ficture_f1_score.isna()
    df["density_label"] = (df.cells_per_mm2 / 1e3).map("{:.1f}".format)

    values = df[metrics].assign(
        MECR_adj=-df.MECR_adj, cells_per_mm2=np.log10(df.cells_per_mm2)
    )
    if ref is None:
        low, high = values[~df.is_control].min(), values[~df.is_control].max()
    df[metrics] = ((values - low) / (high - low)).clip(0, 1)

    df.index = df.index.map(clean_method_name).rename("id")
    return add_scores(df), (fit, low, high)


def build_benchmark_tables(cohorts, group=None, base_path=_constants.BASE_PATH):
    """Build one table per cohort and a combined table averaging their scaled metrics.

    A group table is scaled with the min, max and MECR fit of the all-cells table of the
    same cohort, so its scores are comparable to Fig. 1. Also returns the combined total
    score from pooled scaling, for sensitivity checks.
    """
    tables_used = metric_tables | (group_tables if group else {})
    dfs_by = {
        cohort: {
            table: pd.read_csv(
                Path(base_path) / "metrics" / cohort / folder / f"{table}.csv",
                index_col=0,
            ).query("sample not in ['all', 'all_samples']")
            for table, folder in tables_used.items()
        }
        for cohort in cohorts
    }
    dfs_by["combined"] = {
        t: pd.concat([dfs_by[c][t] for c in cohorts], ignore_index=True)
        for t in tables_used
    }
    tables = {}
    for name, dfs in dfs_by.items():
        tables[name], ref = build_benchmark_table(dfs)
        if group:
            tables[name], _ = build_benchmark_table(dfs, group, ref)
    combined = tables["combined"]
    pooled_score = combined.total_score.copy()
    combined[metrics] = (
        pd.concat([tables[c][metrics] for c in cohorts])
        .groupby(level=0)
        .mean()
        .reindex(combined.index)
    )
    tables["combined"] = add_scores(combined)
    return tables, pooled_score


def plot_benchmark_table(
    df,
    out_file,
    score_label="Overall",
    font_file=Path(_constants.BASE_PATH) / "misc" / "Arial.ttf",
):
    """Plot a benchmark table with funkyheatmap and save it as png and svg."""
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    ro.r["source"](str(Path(__file__).with_suffix(".r")))
    labels = pd.DataFrame(list(_constants.metric_names.items()), columns=["id", "name"])
    with (ro.default_converter + pandas2ri.converter).context():
        ro.r["plot_benchmark_heatmap"](
            df.reset_index(),
            str(out_file),
            labels,
            font_file=str(font_file),
            score_label=score_label,
        )
