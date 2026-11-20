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
    "z": ["low_z_integrity_vs_tissue"],
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
    "z_integrity_all": "z_integrity",
    "z_integrity_tissue": "z_integrity",
}
group_tables = {
    "MECR_score_groups": "marker_gene_metrics",
    "z_integrity_groups": "z_integrity",
}


def z_vs_tissue(z, tissue):
    """Low z integrity share in cells minus the whole-tissue share of the sample (points), mean over samples per method."""
    z = z.query("signal_threshold == 2 and z_integrity_threshold == 0.7")
    diff = z.low_z_integrity_pct - z["sample"].map(
        tissue.set_index("sample").low_z_integrity_pct
    )
    return diff.groupby(z.method).mean().rename("low_z_integrity_vs_tissue")


def add_scores(df):
    """Add total and expression scores as means of the category scores, sorted by total."""
    df["total_score"] = df[list(categories)].mean(axis=1)
    df["expr_score"] = df[["clustering", "marker"]].mean(axis=1)
    return df.sort_values("total_score", ascending=False)


def build_benchmark_table(dfs, group=None, fit=None):
    """Average metrics per method and min-max scale them across non-control methods.

    BAMs are left out of Ficture F1, which has no reliable BAM factor.

    Args:
        dfs: metric tables as loaded by build_benchmark_tables.
        group: key of _constants.celltype_groups to restrict the metrics to, or None for all cells.
            Silhouette is undefined within one cell type and left empty.
        fit: MECR fit on log10 transcripts per cell, fitted on this table if None.

    Returns:
        The scored table (raw metric values in *_raw columns) and its MECR fit.
    """
    general = dfs["general_stats"]
    stats = [
        "volume_final",
        "elongation",
        "intensities_DAPI",
    ]
    all_cells = general.query("cell_type_revised == 'all'").set_index(
        ["sample", "method"]
    )
    marker = dfs["marker_f1_score_cell_type_revised"]
    ficture = dfs["ficture_f1"]
    ficture = ficture[
        ~ficture.cell_type.str.startswith("f1_") & (ficture.cell_type != "BAMs")
    ]
    counts = (
        dfs["assigned_transcript_counts"]
        .groupby(["sample", "method"])[["assigned_count_qced", "total_count"]]
        .sum()
    )
    silhouette = dfs["clustering_score_adata_integrated_cell_type_revised"]
    if group is None:
        cells = all_cells
        transcripts = counts.assigned_count_qced
        mecr = dfs["MECR_score_all"]
        z = dfs["z_integrity_all"]
    else:
        types = _constants.celltype_groups[group]
        rows = general[general.cell_type_revised.isin(types)]
        cells = (
            rows[stats + ["n_counts"]]
            .mul(rows.n_cells, axis=0)
            .assign(n_cells=rows.n_cells, sample=rows["sample"], method=rows.method)
            .groupby(["sample", "method"])
            .sum(min_count=1)
        )
        transcripts = cells.n_counts
        cells = cells.div(cells.n_cells, axis=0).assign(n_cells=cells.n_cells)
        mecr = dfs["MECR_score_groups"].query("group == @group")
        z = dfs["z_integrity_groups"].query("group == @group")
        silhouette = silhouette.iloc[:0]
        marker = marker[
            marker.cell_type.isin(
                {_constants.merged_celltypes.get(t, t) for t in types}
            )
        ]
        ficture = ficture[ficture.cell_type.isin(types)]
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
            z_vs_tissue(z, dfs["z_integrity_tissue"]),
        ],
        axis=1,
    ).drop(index=exclude, errors="ignore")

    df["is_control"] = df.index.str.startswith("Negative_Control")
    if fit is None:
        fit_data = df[~df.is_control].dropna(subset=["MECR", "tpc"])
        fit = np.polyfit(np.log10(fit_data.tpc), fit_data.MECR, 1)
    df["MECR_adj"] = df.MECR - np.polyval(fit, np.log10(df.tpc))

    is_3d = df.index.str.startswith(tuple(_constants.methods_3D))
    is_ficture_seg = df.index.str.startswith("Ficture_segments")
    df.loc[is_ficture_seg, "ficture_f1_score"] = np.nan
    df["na_z"] = is_3d & df.low_z_integrity_vs_tissue.isna()
    df["na_ficture"] = (is_3d | is_ficture_seg) & df.ficture_f1_score.isna()

    df.index = df.index.map(clean_method_name).rename("id")
    df = df.join(df[metrics].add_suffix("_raw"))
    values = df[metrics].assign(
        MECR_adj=-df.MECR_adj,
        low_z_integrity_vs_tissue=-df.low_z_integrity_vs_tissue,
        cells_per_mm2=np.log10(df.cells_per_mm2),
    )
    low, high = values[~df.is_control].min(), values[~df.is_control].max()
    df[metrics] = ((values - low) / (high - low)).clip(0, 1)
    for name, cols in categories.items():
        df[name] = df[cols].mean(axis=1)
    return add_scores(df), fit


def build_benchmark_tables(cohorts, group=None, base_path=_constants.BASE_PATH):
    """Build one table per cohort and a combined table averaging their scaled metrics.

    Group tables use the MECR fit of the all-cells table and are scaled within the group.
    Methods missing in some cohorts are flagged in the "partial" column.
    Also returns the combined total score from pooled scaling, for sensitivity checks.
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
        tables[name], fit = build_benchmark_table(dfs)
        if group:
            tables[name], _ = build_benchmark_table(dfs, group, fit)
    combined = tables["combined"]
    pooled_score = combined.total_score.copy()
    columns = metrics + list(categories) + [f"{m}_raw" for m in metrics]
    combined[columns] = (
        pd.concat([tables[c][columns] for c in cohorts])
        .groupby(level=0)
        .mean()
        .reindex(combined.index)
    )
    tables["combined"] = add_scores(combined)
    indices = [set(tables[c].index) for c in cohorts]
    for table in tables.values():
        table["partial"] = table.index.isin(
            set.union(*indices) - set.intersection(*indices)
        )
    return tables, pooled_score


def plot_benchmark_table(
    df,
    out_file,
    score_label="Overall",
    row_height=1.02,
    width=None,
    height=None,
    note=None,
    font_file=Path(_constants.BASE_PATH) / "misc" / "Arial.ttf",
):
    """Plot a benchmark table with funkyheatmap and save it as png and svg.

    Methods missing in some cohorts are marked with an asterisk. An empty silhouette column is left out.
    width: print width in inches; text is then sized to 6 pt at that width. None keeps the native size.
    height: print height in inches, reached by adding space between rows.
    note: key printed below the table, None for none. With width set, the png is tagged 600 dpi so it places at that width.
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    Path(out_file).parent.mkdir(parents=True, exist_ok=True)
    ro.r["source"](str(Path(__file__).with_suffix(".r")))
    names = _constants.metric_names | {
        "volume_final": "Volume",
    }
    labels = pd.DataFrame(list(names.items()), columns=["id", "name"])
    if df.silhouette_score.isna().all():
        df = df.drop(columns="silhouette_score")
    with (ro.default_converter + pandas2ri.converter).context():
        ro.r["plot_benchmark_heatmap"](
            df.reset_index().assign(
                id=lambda d: d.id + d.partial.map({True: "*", False: ""}),
                density_label=lambda d: (d.cells_per_mm2_raw / 1e3).map(
                    lambda x: f"{x:.1f}" if round(x, 1) < 10 else f"{x:.0f}"
                ),
            ),
            str(out_file),
            labels,
            font_file=str(font_file),
            score_label=score_label,
            row_height=row_height,
            width=ro.NULL if width is None else width,
            height=ro.NULL if height is None else height,
            note=ro.NULL if note is None else note,
        )
    if width is not None:
        from PIL import Image

        png = f"{out_file}.png"
        im = Image.open(png)
        im.load()
        im.save(png, dpi=(600, 600))
