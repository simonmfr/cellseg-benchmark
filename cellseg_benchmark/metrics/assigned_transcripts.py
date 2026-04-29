import pathlib
import functools
import numpy as np
import pandas as pd
import spatialdata as sd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

from .. import _constants

def compute_assigned_transcripts(
    adata,
    *,
    method: str,
    base_path: str | pathlib.Path,
    **kwargs
) -> pd.DataFrame:
    """
    Compute per-gene assigned transcript statistics by combining assigned counts from adata with total detected transcripts from transcript table.

    Args:
        adata: AnnData containing the samples to process; sample labels are expected in `adata.obs["sample"]`.
        method: Segmentation method name, required to access sdata.
        base_path: Base directory, required to access sdata.

    Returns:
        DataFrame with columns
        ["sample", "gene", "assigned_count_qced", "assigned_count_raw",
         "total_count", "pct_assigned_raw", "pct_assigned_qced"].
    """
    df_total = _total_transcripts(adata, base_path=base_path)
    df_assigned = _assigned_transcripts(
        adata=adata,
        method=method,
        base_path=base_path,
    )

    # remove blank genes
    df_total = df_total[~df_total["gene"].str.contains("Blank", case=False, na=False)]
    df_assigned = df_assigned[~df_assigned["gene"].str.contains("Blank", case=False, na=False)]
    assert df_assigned["gene"].nunique() == df_total["gene"].nunique()

    df = df_assigned.merge(df_total, on=["sample", "gene"], how="inner")
    #assert len(df) == len(df_assigned), (
    #    "Mismatch after merging samples from adata and sdata segmentation output. "
    #    "Maybe adata contains samples without segmentation output?"
    #)
    return df.assign(
        pct_assigned_raw=df.assigned_count_raw / df.total_count,
        pct_assigned_qced=df.assigned_count_qced / df.total_count,
    )

def _total_transcripts(
    adata,
    base_path: str | pathlib.Path,
) -> pd.DataFrame:
    """
    Count total detected transcripts per gene per sample from the raw
    transcript table (`points.parquet`) of each sample.

    Returns:
        DataFrame with columns ["sample", "gene", "total_count"].
    """
    base_path = pathlib.Path(base_path)
    samples = sorted(pd.unique(adata.obs["sample"]))
    rows = []

    for sample in samples:
        pfile = (
            base_path
            / "samples"
            / sample
            / "sdata_z3.zarr"
            / "points"
            / f"{sample}_transcripts"
            / "points.parquet"
        )

        if not pfile.exists():
            print(f"points.parquet missing for sample: {sample}")
            continue

        genes = pd.read_parquet(pfile, columns=["gene"])["gene"].astype("string")

        counts = genes.value_counts(dropna=True)

        rows.extend(
            {"sample": sample, "gene": gene, "total_count": int(count)}
            for gene, count in counts.items()
        )

    return pd.DataFrame(rows, columns=["sample", "gene", "total_count"])

def _assigned_transcripts(
    adata,
    method: str,
    base_path: pathlib.Path,
    sample_col: str = "sample",
    layer: str = "counts",
    raw_table_key: str = "table",
) -> pd.DataFrame:
    """
    Compute per-gene assigned transcript counts per sample.

    Uses QC-ed counts from `adata.layers[layer]` and, if available, loads the
    corresponding pre-QC adata from `sdata.zarr`.

    Returns:
        DataFrame with columns
        ["sample", "gene", "assigned_count_qced", "assigned_count_raw"].
    """
    X_qc = adata.layers[layer]
    samples = adata.obs[sample_col].to_numpy()
    genes = adata.var.index.to_numpy()

    records = []
    for sample in pd.unique(samples):
        mask = samples == sample
        qc_counts = X_qc[mask].sum(axis=0)
        qc_counts = qc_counts.A1 if hasattr(qc_counts, "A1") else np.asarray(qc_counts).ravel()

        zarr_path = base_path / "samples" / str(sample) / "results" / method / "sdata.zarr"
        raw_counts = None
        if zarr_path.exists():
            adata_raw = sd.read_zarr(zarr_path)[raw_table_key]
            raw_sum = adata_raw.X.sum(axis=0)
            raw_sum = raw_sum.A1 if hasattr(raw_sum, "A1") else np.asarray(raw_sum).ravel()

            raw_series = pd.Series(raw_sum, index=adata_raw.var.index)
            raw_counts = raw_series.reindex(genes, fill_value=0).to_numpy()
        else:
            print(f"no pre-QC sdata found for sample: {sample}")
        
        if raw_counts is None:
            raw_counts = np.full(len(genes), np.nan)

        records.extend(
            {
                "sample": sample,
                #"method": method,
                "gene": gene,
                "assigned_count_qced": int(qc),
                "assigned_count_raw": (int(raw) if np.isfinite(raw) else np.nan),
            }
            for gene, qc, raw in zip(genes, qc_counts, raw_counts)
        )

    df = pd.DataFrame.from_records(records).astype(
    {"assigned_count_qced": "Int64", "assigned_count_raw": "Int64"}
    )
    
    return df

def plot_assigned_transcripts(cohort: str, boxplot: bool = False, show: bool = True):
    """
    Plot assigned transcript percentages per segmentation method.

    Args:
        cohort: Cohort name used to locate the results CSV.
        boxplot: If True, draw boxplots; otherwise draw bars.
        show: If True, display the figure.
    """
    results_file = (
        pathlib.Path(_constants.BASE_PATH)
        / "metrics"
        / cohort
        / "assigned_transcripts"
        / "assigned_transcript_counts.csv"
    )
    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(results_file, index_col=0)
    df = (
        df.groupby(["sample", "method"], as_index=False)[
            ["assigned_count_qced", "assigned_count_raw", "total_count"]
        ]
        .sum()
    )
    df["pct_assigned_raw"] = df["assigned_count_raw"] / df["total_count"]
    df["pct_assigned_qced"] = df["assigned_count_qced"] / df["total_count"]

    for old, new in _constants.clean_method_names.items():
        df["method"] = df["method"].str.replace(old, new, regex=False)

    # bars: ratio of sums per method
    agg = (
        df.groupby("method", as_index=True)[
            ["assigned_count_raw", "assigned_count_qced", "total_count"]
        ]
        .sum()
    )
    pct_raw = agg["assigned_count_raw"] / agg["total_count"]
    pct_qc  = agg["assigned_count_qced"] / agg["total_count"]

    order    = pct_raw.sort_values().index
    qc_avail = pct_qc.index.intersection(order)

    x = np.arange(len(order))
    w = 0.38

    plt.figure(figsize=(6 + 0.2 * len(order), 4))

    if boxplot:
        raw_data = [df[df["method"] == m]["pct_assigned_raw"].values * 100 for m in order]
        qc_data  = [df[df["method"] == m]["pct_assigned_qced"].values * 100 for m in order]
        bp_kw = dict(widths=w, patch_artist=True, manage_ticks=False,
                     medianprops=dict(color="black"), showfliers=False)
        plt.boxplot(raw_data, positions=x - w / 2,
                    boxprops=dict(facecolor="steelblue"), **bp_kw)
        plt.boxplot(qc_data,  positions=x + w / 2,
                    boxprops=dict(facecolor="lightsteelblue"), **bp_kw)
        h_raw = mpatches.Patch(color="steelblue",      label="All cells (per sample)")
        h_qc  = mpatches.Patch(color="lightsteelblue", label="QCed cells (per sample)")
        plt.xlim(-0.5, len(order) - 0.5)
    else:
        h_raw = plt.bar(x - w / 2, pct_raw[order] * 100, width=w,
                        color="steelblue", label="All cells (pooled)", zorder=1)
        h_qc  = plt.bar(x[order.get_indexer(qc_avail)] + w / 2, pct_qc[qc_avail] * 100,
                        width=w, color="lightsteelblue", label="QCed cells (pooled)", zorder=1)

    idx = pd.Series(x, index=order)

    plt.scatter(
        df["method"].map(idx) - w / 2,
        df["pct_assigned_raw"] * 100,
        s=4,
        color="k",
        alpha=0.2,
        zorder=3,
    )
    plt.scatter(
        df["method"].map(idx) + w / 2,
        df["pct_assigned_qced"] * 100,
        s=4,
        color="k",
        alpha=0.2,
        zorder=3,
    )

    sample_h = plt.Line2D(
        [0], [0],
        marker="o",
        color="black",
        linestyle="none",
        markersize=6,
        alpha=0.6,
    )

    plt.xticks(x, order, rotation=-45, ha="left", va="top")
    plt.ylabel("Assigned Transcripts (%)")
    plt.margins(x=0.02)
    plt.tight_layout()
    plt.legend([h_raw, h_qc, sample_h],
               [h_raw.get_label(), h_qc.get_label(), "Sample"])

    out_file = plot_path / "assigned_transcripts_plot.png"
    plt.savefig(out_file, dpi=300)
    if show:
        plt.show()

def plot_assigned_transcripts_heatmap(
    cohort: str,
    method,
    use_qc: bool = True,
    topn: int = 5,
    show: bool = True,
    save: bool = True,
):
    """
    Plot heatmap of per-gene assigned transcript percentages.

    When selecting multiple methods, genes are ranked by their mean value across all selected methods.

    Args:
        cohort: Cohort.
        method: Segmentation method(s) to plot. Uses original method names (before cleaning).
        use_qc: If True, use QCed percentages; otherwise raw.
        topn: Number of top and bottom genes to include.
        show: If True, display the figure.
        save: If True, save the figure as a PNG in plots directory.

    Returns:
        Path to the saved PNG if save is True, otherwise None.
    """
    results_file = (
        pathlib.Path(_constants.BASE_PATH)
        / "metrics"
        / cohort
        / "assigned_transcripts"
        / "assigned_transcript_counts.csv"
    )
    plot_path = results_file.parent / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(results_file, index_col=0)
    agg = (
        df.groupby(["method", "gene"], as_index=False)[
            ["assigned_count_raw", "assigned_count_qced", "total_count"]
        ].sum()
    )
    agg["pct_raw"] = agg["assigned_count_raw"] / agg["total_count"]
    agg["pct_qc"]  = agg["assigned_count_qced"] / agg["total_count"]

    methods = [method] if isinstance(method, str) else list(method)
    colname = "pct_qc" if use_qc else "pct_raw"

    # build lookup using original method names
    lookup = {}
    display_methods = []
    for m in methods:
        mdf = agg[agg["method"] == m].copy()
        if mdf.empty:
            raise ValueError(f"No rows for method {m!r}")
        lookup[m] = mdf.set_index("gene")[colname]
        display_methods.append(m)

    # rank genes by mean across selected methods, take exactly top/bottom N
    gene_means = pd.Series({
        g: np.nanmean([lookup[m].get(g, np.nan) for m in display_methods])
        for g in agg["gene"].unique()
    }).sort_values()
    candidate_genes = set(gene_means.head(topn).index) | set(gene_means.tail(topn).index)

    # order descending by mean
    ordered_genes = sorted(candidate_genes, key=lambda g: gene_means[g], reverse=True)

    # build matrix
    n_methods = len(display_methods)
    img = np.full((len(ordered_genes), n_methods), np.nan)
    for j, m in enumerate(display_methods):
        for i, g in enumerate(ordered_genes):
            v = lookup[m].get(g, np.nan)
            img[i, j] = np.nan if np.isnan(v) else v * 100

    fig_height = 0.2 * len(ordered_genes) + 0.5
    fig_width  = 0.4 * n_methods + 0.2
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 0.6], wspace=0.1)
    ax  = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])

    cmap = plt.cm.Reds
    norm = plt.Normalize(0, 100)
    im = ax.imshow(img, cmap=cmap, norm=norm, aspect="auto")

    ax.set_yticks(np.arange(len(ordered_genes)))
    ax.set_yticklabels(ordered_genes, fontsize=7)
    ax.set_xticks(np.arange(n_methods))
    ax.set_xticklabels(
        [functools.reduce(lambda n, kv: n.replace(*kv), _constants.clean_method_names.items(), m)
         for m in display_methods],
        fontsize=7, rotation=-45, ha="left", va="top", rotation_mode="anchor",
    )

    for j, m in enumerate(display_methods):
        for i, g in enumerate(ordered_genes):
            v = img[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{int(round(v))}", ha="center", va="center",
                        fontsize=8, color="black" if v < 60 else "white")

    for spine in ("top", "right", "bottom", "left"):
        ax.spines[spine].set_visible(False)

    colorbar = fig.colorbar(im, cax=cax)
    colorbar.outline.set_visible(False)
    colorbar.ax.set_ylabel("Assigned transcripts (%)", fontsize=7)
    colorbar.ax.tick_params(labelsize=7)

    pos = cax.get_position()
    scale = 0.6
    cax.set_position([
        pos.x0,
        pos.y0 + pos.height * (1 - scale) / 2,
        pos.width,
        pos.height * scale,
    ])

    metric = "qc" if use_qc else "raw"
    mtag = (
        functools.reduce(lambda n, kv: n.replace(*kv), _constants.clean_method_names.items(), display_methods[0])
        if len(display_methods) == 1 else f"{n_methods}methods"
    )
    out_file = plot_path / f"gene_pct_assigned_heatmap_{mtag}_{metric}.png"
    if save:
        fig.savefig(out_file, dpi=300, bbox_inches="tight", pad_inches=0.05)
    if show:
        plt.show()
    plt.close(fig)
    return out_file if save else None