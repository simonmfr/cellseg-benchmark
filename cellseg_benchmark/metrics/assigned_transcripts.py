from pathlib import Path
import numpy as np
import pandas as pd
from spatialdata import read_zarr
import matplotlib.pyplot as plt
from matplotlib import gridspec

from cellseg_benchmark import BASE_PATH
from cellseg_benchmark._constants import clean_method_names


def compute_assigned_transcripts(
    adata,
    *,
    method: str,
    base_path: str | Path,
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
    base_path: str | Path,
) -> pd.DataFrame:
    """
    Count total detected transcripts per gene per sample from the raw
    transcript table (`points.parquet`) of each sample.

    Returns:
        DataFrame with columns ["sample", "gene", "total_count"].
    """
    base_path = Path(base_path)
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
    base_path: Path,
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
            adata_raw = read_zarr(zarr_path)[raw_table_key]
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

def plot_assigned_transcripts(cohort: str, show: bool = True):
    results_file = (
        Path(BASE_PATH)
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

    for old, new in clean_method_names.items():
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

    b_raw = plt.bar(
        x - w / 2,
        pct_raw[order] * 100,
        width=w,
        color="steelblue",
        label="All cells (pooled)",
        zorder=1,
    )
    b_qc = plt.bar(
        x[order.get_indexer(qc_avail)] + w / 2,
        pct_qc[qc_avail] * 100,
        width=w,
        color="lightsteelblue",
        label="QCed cells (pooled)",
        zorder=1,
    )

    # scatter: per-sample percentages
    idx = pd.Series(x, index=order)

    plt.scatter(
        df["method"].map(idx) - w / 2,
        df["pct_assigned_raw"] * 100,
        s=7,
        color="k",
        alpha=0.3,
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
    plt.legend(
        [b_raw, b_qc, sample_h],
        ["All cells (pooled)", "QCed cells (pooled)", "Sample"],
    )

    out_file = plot_path / "assigned_transcripts_plot.png"
    plt.savefig(out_file, dpi=300)
    if show:
        plt.show()

def plot_assigned_transcripts_heatmap(cohort: str,
                                      method,
                                      use_qc: bool = True,
                                      topn: int = 5,
                                      show: bool = True,
                                      save: bool = True):
    """
    Plot heatmap of per-gene assigned transcript percentages.
    
    Args:
        cohort: Cohort.
        method: Segmentation method(s) to plot.
        use_qc: If True, use QCed percentages; otherwise raw.
        topn: Number of top and bottom genes to include.
        show: If True, display the figure.
        save: If True, save the figure as a PNG in plots directory.
    
    Returns:
        Path to the saved PNG if save is True, otherwise None.
    """
    results_file = (
        Path(BASE_PATH)
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

    for old, new in clean_method_names.items():
        agg["method"] = agg["method"].str.replace(old, new, regex=False)

    methods = [method] if isinstance(method, str) else list(method)
    colname = "pct_qc" if use_qc else "pct_raw"

    vals_list, genes_list = [], []
    display_methods = []

    for m in methods:
        m_display = clean_method_names.get(m, m)
        mdf = agg[agg["method"] == m_display].copy()
        if mdf.empty:
            raise ValueError(f"No rows for method {m_display!r}")
        mdf = mdf.sort_values(colname)
        sel = pd.concat([mdf.head(topn), mdf.tail(topn)]).sort_values(
            colname, ascending=False
        )
        genes_list.append(sel["gene"].to_numpy())
        vals_list.append(sel[colname].to_numpy() * 100)
        display_methods.append(m_display)

    n_methods = len(display_methods)
    max_genes = max(len(g) for g in genes_list)

    img = np.full((max_genes, n_methods), np.nan)
    for j, vals in enumerate(vals_list):
        img[: len(vals), j] = vals

    fig_height = 0.2 * max_genes + 0.5
    fig_width  = 0.4 * n_methods + 0.2
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 0.6], wspace=0.1)
    ax  = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])

    cmap = plt.cm.Reds
    norm = plt.Normalize(0, 100)
    im = ax.imshow(img, cmap=cmap, norm=norm, aspect="auto")

    ax.set_yticks(np.arange(len(genes_list[0])))
    ax.set_yticklabels(genes_list[0], fontsize=7)
    ax.set_xticks(np.arange(n_methods))
    ax.set_xticklabels(display_methods, fontsize=7, rotation=-45,
                       ha="left", va="top", rotation_mode="anchor")

    for j, vals in enumerate(vals_list):
        for i, v in enumerate(vals):
            ax.text(
                j, i, f"{int(round(v))}",
                ha="center", va="center",
                fontsize=8,
                color="black" if v < 60 else "white",
            )

    for spine in ("top", "right", "bottom", "left"):
        ax.spines[spine].set_visible(False)

    cb = fig.colorbar(im, cax=cax)
    cb.outline.set_visible(False)
    cb.ax.set_ylabel("Assigned transcripts (%)", fontsize=7)
    cb.ax.tick_params(labelsize=7)

    pos = cax.get_position()
    scale = 0.6
    cax.set_position([
        pos.x0,
        pos.y0 + pos.height * (1 - scale) / 2,
        pos.width,
        pos.height * scale,
    ])

    metric = "qc" if use_qc else "raw"
    if len(display_methods) == 1:
        mtag = display_methods[0]
    else:
        mtag = f"{len(display_methods)}methods"
    out_file = plot_path / f"gene_pct_assigned_heatmap_{mtag}_{metric}.png"
    if save:
        fig.savefig(
            out_file,
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.05,
        )
    if show:
        plt.show()
    plt.close(fig)
    return out_file if save else None