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


def add_scores(df):
    """Add total and expression scores as means of category means, sorted by total."""
    category_scores = pd.DataFrame(
        {name: df[cols].mean(axis=1) for name, cols in categories.items()}
    )
    df["total_score"] = category_scores.mean(axis=1)
    df["expr_score"] = category_scores[["clustering", "marker"]].mean(axis=1)
    return df.sort_values("total_score", ascending=False)


def build_benchmark_table(dfs):
    """Average metrics per method and min-max scale them across non-control methods."""
    general = dfs["general_stats"].query("cell_type_revised == 'all'")
    counts = (
        dfs["assigned_transcript_counts"]
        .groupby(["sample", "method"])[["assigned_count_qced", "total_count"]]
        .sum()
    )
    counts["pct_assigned_qced"] = counts.assigned_count_qced / counts.total_count
    counts["tpc"] = (
        counts.assigned_count_qced / general.set_index(["sample", "method"]).n_cells
    )
    density = dfs["cell_density"].groupby("method")[["n_cells", "tissue_mm2"]].sum()
    marker_f1 = (
        dfs["marker_f1_score_cell_type_revised"]
        .groupby(["sample", "method", "cell_type"])
        .f1_score.mean()
        .groupby(["sample", "method"])
        .mean()
        .groupby("method")
        .mean()
    )

    df = pd.concat(
        [
            general.groupby("method")[
                [
                    "volume_final",
                    "elongation",
                    "intensities_DAPI",
                    "Ovrlpy_stats_mean_integrity",
                ]
            ].mean(),
            counts.groupby("method")[["pct_assigned_qced", "tpc"]].mean(),
            dfs["clustering_score_adata_integrated_cell_type_revised"]
            .groupby("method")[["silhouette_score"]]
            .mean(),
            dfs["MECR_score_all"]
            .query("gene1 == 'all'")
            .groupby("method")[["MECR"]]
            .mean(),
            marker_f1.rename("marker_f1"),
            dfs["ficture_f1"]
            .query("cell_type == 'f1_macro_all'")
            .groupby("method")
            .f1.mean()
            .rename("ficture_f1_score"),
            (density.n_cells / density.tissue_mm2).rename("cells_per_mm2"),
        ],
        axis=1,
    ).drop(index=exclude, errors="ignore")

    df["is_control"] = df.index.str.startswith("Negative_Control")
    fit_data = df[~df.is_control].dropna(subset=["MECR", "tpc"])
    fit = np.polyfit(np.log10(fit_data.tpc), fit_data.MECR, 1)
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
    low, high = values[~df.is_control].min(), values[~df.is_control].max()
    df[metrics] = ((values - low) / (high - low)).clip(0, 1)

    df.index = df.index.map(clean_method_name).rename("id")
    return add_scores(df)


def build_benchmark_tables(cohorts, base_path=_constants.BASE_PATH):
    """Build one table per cohort and a combined table averaging their scaled metrics.

    Also returns the combined total score from pooled scaling, for sensitivity checks.
    """
    dfs_by = {
        cohort: {
            table: pd.read_csv(
                Path(base_path) / "metrics" / cohort / folder / f"{table}.csv",
                index_col=0,
            ).query("sample not in ['all', 'all_samples']")
            for table, folder in metric_tables.items()
        }
        for cohort in cohorts
    }
    dfs_by["combined"] = {
        t: pd.concat([dfs_by[c][t] for c in cohorts], ignore_index=True)
        for t in metric_tables
    }
    tables = {name: build_benchmark_table(dfs) for name, dfs in dfs_by.items()}
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
