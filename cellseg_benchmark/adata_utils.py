import logging
import math
from os.path import isfile, join
from typing import List, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from anndata import AnnData, concat
from matplotlib.path import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable
from spatialdata import SpatialData
from tqdm import tqdm

from cellseg_benchmark.cell_annotation_utils import cell_type_colors


def merge_adatas(
    sdatas: List[Tuple[str, SpatialData]],
    seg_method: str,
    sample_paths_file: dict,
    logger: logging.Logger = None,
    plot_qc_stats=False,
    save_path=None,
) -> AnnData:
    """Merge AnnData objects extracted from multiple SpatialData objects into a single AnnData.

    Each SpatialData object contains results of multiple segmentation approaches (for a single
    respective sample). This function extracts AnnData objects for a given segmentation method key,
    annotates them with metadata, optionally exports QC visualizations, and merges them into a
    single AnnData object. No integration or batch effect correction is performed.

    Args:
        sdatas (List[Tuple[str, SpatialData]]): List of tuples, each containing a unique
            sample name (str) and a SpatialData object. Each SpatialData is expected to
            have an AnnData accessible via `sdata[f"adata_{seg_method}"]`.
        seg_method (str): Segmentation method key.
        sample_paths_file (dict): Dict that maps sample_name to Merscope output data path.
        logger (logging.Logger, optional): Logger instance for informational and warning
            messages. If None, logging is disabled. Defaults to None.
        plot_qc_stats (bool, optional): If True, triggers QC plotting via the `plot_qc` function internally. Defaults to False.
        save_path (str, optional): File path or directory where QC plots will be saved
            if `plot_qc_stats` is True. Defaults to None.

    Returns:
        AnnData: A merged AnnData object containing concatenated cells from all input datasets.
    """
    adatas = []

    y_limits = [0, 0, 0, 0]
    if logger:
        logger.info(f"Merging adatas of {seg_method}")
    for name, sdata in tqdm(sdatas):
        if f"adata_{seg_method}" not in sdata.tables.keys():
            if logger:
                logger.warning(f"Skipping {name}. No such key: {seg_method}")
            continue
        adata = sdata[f"adata_{seg_method}"]
        samples = name.split("_")
        adata.obs["cohort"] = samples[0]
        adata.obs["slide"] = samples[1]
        adata.obs["region"] = samples[2]
        adata.obs["condition"] = sample_paths_file[name].split("/")[-2].split("-")[-2]
        if isinstance(adata.X, np.ndarray):
            adata.X = sp.csr_matrix(adata.X, dtype=np.float32)
        adata.obs["n_counts"] = adata.X.sum(axis=1)
        adata.obs["n_genes"] = adata.X.count_nonzero(axis=1)
        adata.obs["sample"] = name
        adatas.append(adata)

        if plot_qc_stats:
            y_limits[0] = max(
                y_limits[0], max(np.histogram(adata.obs["n_counts"], bins=60)[0])
            )
            y_limits[1] = max(
                y_limits[1],
                max(
                    np.histogram(
                        adata.obs["n_counts"][adata.obs["n_counts"] < 1000], bins=60
                    )[0]
                ),
            )
            y_limits[2] = max(
                y_limits[2], max(np.histogram(adata.obs["n_genes"], bins=60)[0])
            )
            y_limits[3] = max(
                y_limits[3], max(np.histogram(adata.obs["volume"], bins=100)[0])
            )

    adata = concat(adatas, join="outer", merge="first")
    for i in set(["cell_id", "cell_id_x", "cell_id_y"]) & set(adata.obs.columns):
        del adata.obs[i]
    adata.obs["cell_type_mmc_raw_revised"] = pd.Categorical(
        adata.obs["cell_type_mmc_raw_revised"], categories=list(cell_type_colors.keys())
    )
    adata.obs["cell_type_mmc_raw_revised"] = adata.obs[
        "cell_type_mmc_raw_revised"
    ].cat.remove_unused_categories()

    # set colors in adata.uns
    colors = [
        cell_type_colors[cat]
        for cat in adata.obs["cell_type_mmc_raw_revised"].cat.categories
    ]
    adata.uns["cell_type_mmc_raw_revised_colors"] = colors
    adata.obs_names_make_unique()
    if logger:
        n = len(adata.obs["sample"].unique())
        logger.info(f"{seg_method}: # of cells: {len(adata)}, # of samples: {n}")
    if plot_qc_stats:
        plot_qc(adata, save_path, y_limits, logger)
    return adata


def plot_qc(adata: AnnData, save_dir, y_limits, logger) -> None:
    """Plot general quality control stats.

    Args:
        adata: adata object
        save_dir: Plotting directory
        y_limits: Limits for plotting
        logger: logger object

    Returns:
        None
    """
    if logger:
        logger.info("Plotting QC results")
    sample_names = adata.obs["sample"].unique()

    # General Stats
    fig, axs = plt.subplots(
        len(sample_names),
        4,
        figsize=(16, len(sample_names) * 4),
        gridspec_kw={"wspace": 0.4},
        squeeze=False,
    )
    for i, name in enumerate(sample_names):
        adata_tmp = adata[adata.obs["sample"] == name]

        for ax, data, bins, y_limit in zip(
            axs[i, :],
            [
                adata_tmp.obs["n_counts"],
                adata_tmp.obs["n_counts"][adata_tmp.obs["n_counts"] < 1000],
                adata_tmp.obs["n_genes"],
                adata_tmp.obs["volume"],
            ],
            [60, 60, 60, 100],
            y_limits,
        ):
            sns.histplot(data, kde=False, bins=bins, ax=ax, zorder=2, alpha=1)
            ax.grid(True, zorder=1)
            ax.set_ylim(0, y_limit)  # Set consistent y-axis limit for each column
            for patch in ax.patches:
                patch.set_edgecolor("none")  # Remove border lines

    pad = 5
    cols = [
        "Number of Counts per Cell",
        "Number of Counts per Cell (until n=1000)",
        "Number of Genes per Cell",
        "Volume",
    ]
    rows = sample_names

    for ax, col in zip(axs[0, :], cols):
        ax.annotate(
            col,
            xy=(0.5, 1),
            xytext=(0, pad),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="center",
            va="baseline",
        )

    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad - pad, 0),
            xycoords=ax.yaxis.label,
            textcoords="offset points",
            size="large",
            ha="right",
            va="center",
            rotation="vertical",
        )

    fig.tight_layout()
    fig.savefig(join(save_dir, "qc_general_stats_cells.png"))
    plt.close()

    # Cell Qualities
    fig, axs = plt.subplots(
        len(sample_names),
        1,
        figsize=(9, len(sample_names) * 8),
        gridspec_kw={"wspace": 0.4},
    )
    for ax, name in zip(axs, sample_names):
        adata_tmp = adata[adata.obs["sample"] == name]
        sc.pl.scatter(
            adata_tmp,
            x="n_counts",
            y="n_genes",
            color="volume",
            ax=ax,
            show=False,
            legend_loc="none",
        )

        default_cax = fig.axes[-1]
        if default_cax is not ax:
            fig.delaxes(default_cax)

        im = ax.collections[-1]
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label("volume")

    for ax, row in zip(axs, rows):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad - pad, 0),
            xycoords=ax.yaxis.label,
            textcoords="offset points",
            size="large",
            ha="right",
            va="center",
            rotation="vertical",
        )

    fig.tight_layout()
    fig.savefig(join(save_dir, "qc_scatterplot_cells.png"))
    plt.close()

    # General Stats Genes
    fig, axs = plt.subplots(
        len(sample_names),
        3,
        figsize=(18, len(sample_names) * 5),
        gridspec_kw={"wspace": 0.4},
    )
    for i, name in enumerate(sample_names):
        adata_tmp = adata[adata.obs["sample"] == name]
        adata_tmp.var["n_counts"] = adata_tmp.X.sum(axis=0).T
        log10_counts = np.log10(adata_tmp.var["n_counts"] + 1)

        sns.violinplot(y=log10_counts, ax=axs[i, 0], zorder=2)
        axs[i, 0].set_ylabel("log10(n_counts+1)")
        axs[i, 0].grid(True, alpha=0.5, zorder=0)

        sns.histplot(log10_counts, kde=True, bins=25, ax=axs[i, 1], zorder=2)
        axs[i, 1].set_xlabel("log10(n_counts+1)")
        axs[i, 1].grid(True, alpha=0.5, zorder=0)

        sns.ecdfplot(adata_tmp.var["n_counts"], log_scale=True, ax=axs[i, 2], zorder=2)
        axs[i, 2].grid(True, alpha=0.5, zorder=0)

    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad - pad, 0),
            xycoords=ax.yaxis.label,
            textcoords="offset points",
            size="large",
            ha="right",
            va="center",
            rotation="vertical",
        )

    fig.tight_layout()
    fig.savefig(join(save_dir, "qc_general_stats_genes.png"))
    plt.close()

    # Highly expressed genes
    fig, axs = plt.subplots(
        len(sample_names),
        1,
        figsize=(6, len(sample_names) * 4),
    )
    for ax, name in zip(axs, sample_names):
        adata_tmp = adata[adata.obs["sample"] == name]
        sc.pl.highest_expr_genes(adata_tmp, n_top=20, ax=ax, show=False)

    for ax, row in zip(axs, rows):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad - pad, 0),
            xycoords=ax.yaxis.label,
            textcoords="offset points",
            size="large",
            ha="right",
            va="center",
            rotation="vertical",
        )

    fig.tight_layout()
    fig.savefig(join(save_dir, "qc_highly_expressed_genes.png"))
    plt.close()


def filter_spatial_outlier_cells(
    adata: AnnData,
    data_dir: str,
    sample_paths_file: dict,
    save_path: str,
    logger: logging.Logger = None,
    remove_outliers: bool = True,
) -> AnnData:
    """Remove cells outside main tissue using manually defined spatial outlier regions.

    Reads outlier coordinates from cell_outlier_coordinates.csv files (generated by
    lasso selection in 10X Explorer) and marks cells within these polygons as spatial outliers.
    Exports qc_spatial_outlier_cells.png showing outlier cells in red.

    Args:
        adata (AnnData): AnnData object with spatial coordinates in micron units
            in obsm["spatial"].
        data_dir (str): Directory containing sample subdirectories with CSV files.
        sample_paths_file (dict): Dict that maps sample_name to Merscope output data path.
        save_path (str): Directory for saving plots.
        logger (logging.Logger, optional): Optional logger object. Defaults to None.
        remove_outliers (bool, optional): If True, remove outlier cells from adata.
            If False, only mark outliers without filtering. Defaults to True.

    Returns:
        AnnData: Modified adata with "spatial_outlier" column added to obs.
            If remove_outliers=True, outlier cells are removed.

    Note:
        Requires cell_outlier_coordinates.csv with 4 header lines starting with #
        and columns: Selection, X, Y
    """
    n_cols = 3
    samples = adata.obs["sample"].unique()
    n_samples = len(samples)
    n_rows = math.ceil(n_samples / n_cols)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4 * n_cols, 4 * n_rows),
        squeeze=False,
        sharex=False,
        sharey=False,
    )

    adata.obs["spatial_outlier"] = False

    for idx, sample in enumerate(samples):
        ax = axes[idx // n_cols][idx % n_cols]

        mask = adata.obs["sample"] == sample
        coords = adata.obsm["spatial"][mask]

        # Subsample for plotting if needed
        sample_data = adata[mask]
        if sample_data.n_obs > int(1.5e5):
            sample_data = sc.pp.subsample(
                sample_data, n_obs=int(1.5e5), random_state=42, copy=True
            )
        plotting_coords = sample_data.obsm["spatial"]

        csv_path = join(data_dir, "samples", sample, "cell_outlier_coordinates.csv")
        if not isfile(csv_path):
            if logger:
                logger.info(
                    f"[{sample}] Missing polygon file: {csv_path}. Draw outliers in 10X Explorer if necessary."
                )
            adata.obs.loc[mask, "spatial_outlier"] = False

            # Plot all cells as non-outliers
            ax.scatter(
                plotting_coords[:, 0],
                plotting_coords[:, 1],
                c="lightgrey",
                s=max(0.3, min(0.7, 30000 / len(plotting_coords))),
                alpha=0.75,
                edgecolors="none",
            )
            continue

        coords_df = pd.read_csv(csv_path, skiprows=4)  # skip comment header lines
        if list(coords_df.columns[:3]) != ["Selection", "X", "Y"]:
            raise ValueError(
                f"[{sample}] Unexpected columns in cell_outlier_coordinates.csv: {list(coords_df.columns[:3])}; expected ['Selection', 'X', 'Y']"
            )

        transform = np.loadtxt(
            join(
                sample_paths_file[sample],
                "images",
                "micron_to_mosaic_pixel_transform.csv",
            )
        ).reshape(3, 3)
        coords_df[["X", "Y"]] -= transform[[0, 1], 2] / transform[[0, 1], [0, 1]]

        # Detect outliers on full dataset
        outliers = np.zeros(len(coords), dtype=bool)
        for sel in coords_df["Selection"].unique():
            poly = coords_df[coords_df["Selection"] == sel][["X", "Y"]].values
            if len(poly) > 2:
                outliers |= Path(poly).contains_points(coords)

        adata.obs.loc[mask, "spatial_outlier"] = outliers

        # Detect outliers on plotting dataset for visualization
        plotting_outliers = np.zeros(len(plotting_coords), dtype=bool)
        for sel in coords_df["Selection"].unique():
            poly = coords_df[coords_df["Selection"] == sel][["X", "Y"]].values
            if len(poly) > 2:
                plotting_outliers |= Path(poly).contains_points(plotting_coords)

        colors = pd.Categorical(
            np.where(plotting_outliers, "red", "lightgrey"),
            categories=["lightgrey", "red"],
        )

        ax.scatter(
            plotting_coords[:, 0],
            plotting_coords[:, 1],
            c=colors,
            s=max(0.3, min(0.7, 30000 / len(plotting_coords))),
            alpha=0.75,
            edgecolors="none",
        )

        # Draw polygons
        # for sel in coords_df["Selection"].unique():
        #     poly = coords_df[coords_df["Selection"] == sel][["X", "Y"]].values
        #     if len(poly) > 2:
        #         ax.add_patch(Polygon(poly, fill=False, edgecolor="red", linewidth=0.8))

        ax.set_title(sample, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    for ax in axes.flat[n_samples:]:
        ax.set_visible(False)

    fig.suptitle("Spatial outlier cells by sample", fontsize=14, y=1)
    fig.tight_layout()
    fig.savefig(
        join(save_path, "qc_spatial_outlier_cells.png"), dpi=200, bbox_inches="tight"
    )
    plt.close()

    if logger:
        total_before = adata.n_obs
        n_spatial_outlier = adata.obs["spatial_outlier"].sum()

        logger.info(f"# total cells before filtering spatial outliers: {total_before}")
        logger.info(f"# spatial_outlier:       {n_spatial_outlier}")

    if remove_outliers:
        adata = adata[~adata.obs["spatial_outlier"]].copy()

    if logger:
        total_after = adata.n_obs
        logger.info(f"# total cells after filtering spatial outliers:  {total_after}")

    return adata


def filter_low_quality_cells(
    adata: AnnData,
    save_path: str,
    min_counts: int = 25,
    min_genes: int = 5,
    logger: logging.Logger = None,
    remove_outliers: bool = True,
) -> AnnData:
    """Flags two types of problematic cells in `adata.obs`.

    - 'low_quality_cell': cells with fewer than `min_counts` counts or `min_genes` genes.
    - 'volume_outlier_cell': cells with volume < 100 or > 3× median volume.

    If `remove_outliers` is True:
    - Volume outliers are removed.
    - `scanpy`-based filtering is applied to remove low-quality cells.

    Plots of flagged cells are saved in `save_path`.

    Args:
        adata (AnnData): AnnData object with spatial coordinates in micron units
            in obsm["spatial"].
        save_path (str): Directory for saving plots.
        min_counts (int): Minimum number of counts per cell to retain.
        min_genes (int): Minimum number of genes per cell to retain.
        logger (logging.Logger, optional): Optional logger object. Defaults to None.
        remove_outliers (bool, optional): If True, remove outlier cells from adata.
            If False, only mark outliers without filtering. Defaults to True.

    Returns:
        AnnData: Modified adata with "spatial_outlier" column added to obs.
            If remove_outliers=True, outlier cells are removed.
    """
    adata.obs["low_quality_cell"] = (adata.obs["n_counts"] < min_counts) | (
        adata.obs["n_genes"] < min_genes
    )

    metric, n = "volume", 3
    adata.obs["volume_outlier_cell"] = (
        adata.obs[metric] > n * np.median(adata.obs[metric])
    ) | (adata.obs[metric] < 100)

    def _plot_flag(flag, fname):
        n_cols = 3
        samples = adata.obs["sample"].unique()
        n_samples = len(samples)
        n_rows = math.ceil(n_samples / n_cols)

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(4 * n_cols, 4 * n_rows),
            squeeze=False,
            sharex=False,
            sharey=False,
        )

        for idx, sample in enumerate(samples):
            ax = axes[idx // n_cols][idx % n_cols]

            mask = adata.obs["sample"] == sample
            sample_data = adata[mask]

            if sample_data.n_obs > int(1.5e5):
                sample_data = sc.pp.subsample(
                    sample_data, n_obs=int(1.5e5), random_state=42, copy=True
                )

            coords = sample_data.obsm["spatial"]

            outliers = sample_data.obs[flag].values
            colors = pd.Categorical(
                np.where(outliers, "red", "lightgrey"), categories=["lightgrey", "red"]
            )

            ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=colors,
                s=max(0.3, min(0.7, 30000 / len(coords))),
                alpha=0.75,
                edgecolors="none",
            )

            ax.set_title(sample, fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)

        # Hide unused subplots
        for ax in axes.flat[n_samples:]:
            ax.set_visible(False)

        fig.suptitle(f"{flag.replace('_', ' ').title()} by sample", fontsize=14, y=1)
        fig.tight_layout()
        fig.savefig(join(save_path, fname), dpi=200, bbox_inches="tight")
        plt.close(fig)

    _plot_flag("low_quality_cell", "qc_low_quality_cells.png")
    _plot_flag("volume_outlier_cell", "qc_volume_outlier_cells.png")

    if logger:
        total_before = adata.n_obs
        n_low_quality = adata.obs["low_quality_cell"].sum()
        n_volume_outlier = adata.obs["volume_outlier_cell"].sum()

        logger.info(
            f"# total cells before filtering low-quality or volume outliers: {total_before}"
        )
        logger.info(f"# low_quality_cell:             {n_low_quality}")
        logger.info(f"# volume_outlier_cell:       {n_volume_outlier}")

    if remove_outliers:
        adata = adata[~adata.obs["volume_outlier_cell"]].copy()
        sc.pp.filter_cells(adata, min_counts=min_counts)
        sc.pp.filter_cells(adata, min_genes=min_genes)

    if logger:
        total_after = adata.n_obs
        logger.info(
            f"# total cells after filtering low-quality or volume outliers:  {total_after}"
        )

    return adata


def filter_genes(adata, save_path, logger=None):
    """Filter adata after genes.

    Args:
        adata: AnnData
        save_path: save path for plots
        logger: logger object

    Returns:
        adata with filtered genes
    """
    if logger:
        logger.info(f"# genes before filtering: {adata.n_vars}")
    sc.pp.filter_genes(adata, min_cells=10)
    if logger:
        logger.info(f"# genes after filtering: {adata.n_vars}")

    fig, axs = plt.subplots(1, 4, figsize=(18, 4), gridspec_kw={"wspace": 0.4})
    for ax, data, bins in zip(
        axs,
        [
            adata.obs["n_counts"],
            adata.obs["n_counts"][adata.obs["n_counts"] < 1000],
            adata.obs["n_genes"],
            adata.obs["volume"],
        ],
        [100, 100, 100, 100],  # bins
    ):
        sns.histplot(data, kde=False, bins=bins, ax=ax, zorder=2, alpha=1)
        ax.grid(True, zorder=1)
        for patch in ax.patches:
            patch.set_edgecolor("none")
    fig.savefig(join(save_path, "qc_general_stats.png"))
    plt.close()
    return adata


def normalize_counts(
    adata: AnnData,
    save_path: str,
    seg_method: str,
    *,
    target_sum: int = 250,
    logger=None,
) -> AnnData:
    """Volume‑based normalisation of raw UMI counts (as in Allen et al. Cell, 2023).

    Layers added:
    - counts: raw counts
    - volume_log1p_norm: log1p‐transformed volume‑normalised counts
    - zscore: z‑scored volume_log1p_norm
    - librarysize_log1p_norm: log1p library‑size normalised counts (comparison only)

    Args:
        adata (AnnData): AnnData with raw integer counts in ``adata.X`` and a
            ``volume`` column in ``adata.obs``.
        save_path (str): Directory for diagnostic plots.
        seg_method (str): Name of segmentation method for processing logic.
        target_sum (int, optional): Target sum per cell after rescaling.
            Defaults to 250 as in Allen et al.
        logger (logging.Logger, optional): Python logging instance. Defaults to None.

    Returns:
        AnnData: ``adata`` with the layers above and outlier cells removed.
    """
    if logger:
        logger.info("Normalizing counts...")

    if sp.issparse(adata.X) and not isinstance(adata.X, sp.csr_matrix):
        adata.X = adata.X.tocsr()

    if (
        not np.issubdtype(adata.X.dtype, np.integer)
        and "proseg" not in seg_method.lower()
    ):  # exception for proseg: counts are non-integer posterior expectations
        raise TypeError(
            f"adata.X must contain integer counts, found instead: {adata.X.dtype}"
        )

    # 1 Volume normalisation
    adata.layers["counts"] = adata.X
    inv_vol = (1.0 / adata.obs["volume"].to_numpy()).astype("float32")
    adata.layers["volume_norm"] = sp.csr_matrix(adata.X, dtype=np.float32).multiply(
        inv_vol[:, None]
    )

    # 2 Remove outlier cells (1–99 %ile)
    row_sums = np.ravel(adata.layers["volume_norm"].sum(1))
    lo, hi = np.percentile(row_sums, [1, 99])
    mask = (row_sums > lo) & (row_sums < hi)
    if logger:
        logger.info(
            f"Cells before/after outlier removal during normalization: {adata.n_obs} → {mask.sum()}"
        )

    for layer_name in adata.layers.keys():
        if sp.issparse(adata.layers[layer_name]) and not isinstance(
            adata.layers[layer_name], sp.csr_matrix
        ):
            adata.layers[layer_name] = adata.layers[layer_name].tocsr()
    adata = adata[mask].copy()

    # diagnostic plot
    sns.histplot(row_sums, kde=False, edgecolor="none", color="steelblue")
    for p in (lo, hi):
        plt.axvline(p, c="red", lw=1)
    plt.savefig(join(save_path, "normalize_count_outlier_cells.png"))
    plt.close()

    # 3 Rescale per cell to ``target_sum``
    row_sums = np.ravel(adata.layers["volume_norm"].sum(1))
    adata.layers["volume_norm"] = (
        adata.layers["volume_norm"]
        .multiply((target_sum / row_sums).astype("float32")[:, None])
        .astype("float32")
    )

    # 4 Log‑transform and z‑score
    adata.layers["volume_norm"] = sp.csr_matrix(
        adata.layers["volume_norm"], dtype=np.float32
    )
    adata.layers[("volume_log1p_norm")] = sc.pp.log1p(
        adata.layers["volume_norm"], copy=True
    )
    adata.layers["zscore"] = sc.pp.scale(
        adata.layers["volume_log1p_norm"], zero_center=True, max_value=None
    )

    # 5 Library‑size normalisation (comparison only)
    libnorm = sc.pp.normalize_total(
        adata, layer="counts", target_sum=None, inplace=False
    )["X"]
    libnorm = sp.csr_matrix(libnorm, dtype=np.float32)
    adata.layers["librarysize_log1p_norm"] = sc.pp.log1p(libnorm, copy=True)

    # 6 Compare distributions
    fig, axes = plt.subplots(1, 4, figsize=(18, 4))
    layers = [
        ("counts", "Raw counts"),
        ("volume_log1p_norm", "Volume log1p-normalized"),
        ("zscore", "Volume z‑scored"),
        ("librarysize_log1p_norm", "Lib‑size log1p-normalized"),
    ]
    for ax, (layer, title) in zip(axes, layers):
        sns.histplot(
            adata.layers[layer].sum(1),
            bins=100,
            ax=ax,
            edgecolor="none",
            color="steelblue",
            alpha=1,
            legend=False,
            zorder=2,
        )
        ax.set_title(title)
        ax.grid(True, zorder=-1)
    max_ylim = max(ax.get_ylim()[1] for ax in axes)
    for ax in axes:
        ax.set_ylim(0, max_ylim)
    plt.tight_layout()
    plt.savefig(join(save_path, "normalize_count_distributions.png"))
    plt.close()

    assert round(adata.layers["volume_norm"].sum(1).mean()) == target_sum

    return adata


def dimensionality_reduction(adata: AnnData, save_path: str, logger=None) -> AnnData:
    """Run PCA and multiple UMAP projections on the z-scored layer of the input AnnData object.

    Exports PCA.png and comparing_UMAPs.png.

    Args:
        adata (AnnData): Annotated data matrix with a 'zscore' layer and relevant metadata in .obs.
        save_path (str): Directory path to save PCA and UMAP plots.
        logger (optional): Logger instance for status messages.

    Returns:
        AnnData
            Updated AnnData object with PCA and multiple UMAP embeddings added to .obsm.
    """
    adata.X = adata.layers["zscore"]
    sc.settings.figdir = save_path

    if logger:
        logger.info("Dimensionality reduction: PCA")
    sc.pp.pca(adata, svd_solver="arpack")
    fig, axs = plt.subplots(1, 2, figsize=(20, 8), gridspec_kw={"wspace": 0.25})
    with plt.ioff():
        sc.pl.pca_scatter(adata, color="n_counts", ax=axs[0], show=False)
    var_ratio = adata.uns["pca"]["variance_ratio"]

    axs[1].plot(np.arange(1, len(var_ratio) + 1), var_ratio, marker="o", linestyle="-")
    axs[1].set_yscale("log")
    axs[1].set_xlabel("Principal component")
    axs[1].set_ylabel("Variance ratio")
    axs[1].set_title("PCA variance ratio")
    fig.savefig(join(save_path, "PCA.png"))
    plt.close(fig)

    pt_size_umap = 320000 / adata.shape[0]

    fig, axs = plt.subplots(
        3, 3, figsize=(21, 19), gridspec_kw={"wspace": 0.6, "hspace": 0.3}
    )

    color_schemes = ["sample", "condition", "cell_type_mmc_raw_revised"]
    umap_configs = [
        {"n_neighbors": 10, "n_pcs": 20, "key": "X_umap_10_20"},
        {"n_neighbors": 15, "n_pcs": 40, "key": "X_umap_15_40"},
        {"n_neighbors": 20, "n_pcs": 50, "key": "X_umap_20_50"},
    ]

    for row, config in enumerate(umap_configs):
        if logger:
            logger.info(
                f"Dimensionality reduction: UMAP with n_neighbors={config['n_neighbors']} and n_pcs={config['n_pcs']}"
            )

        sc.pp.neighbors(adata, n_neighbors=config["n_neighbors"], n_pcs=config["n_pcs"])
        sc.tl.umap(
            adata, neighbors_key="neighbors", key_added=config["key"], n_components=2
        )

        adata.uns.pop("neighbors", None)
        adata.obsp.pop("distances", None)
        adata.obsp.pop("connectivities", None)

        for col, color_scheme in enumerate(color_schemes):
            with plt.ioff():
                with plt.style.context("default"):
                    with mpl.rc_context({"figure.figsize": (7, 7)}):
                        sc.pl.embedding(
                            adata,
                            ax=axs[row, col],
                            basis=config["key"],
                            size=pt_size_umap,
                            color=color_scheme,
                            title=f"UMAP unintegrated (n_neigh={config['n_neighbors']}, n_pcs={config['n_pcs']}) \n {color_scheme}",
                            show=False,
                        )

    plt.savefig(
        join(save_path, "UMAP_unintegrated_comparison.png"),
        dpi=200,
        bbox_inches="tight",
    )
    plt.close()

    return adata


def integration_harmony(
    adata: AnnData,
    batch_key: str,
    save_path: str,
    logger: logging.Logger = None,
    n_neighbors: int = 20,
    n_pcs: int = 50,
) -> AnnData:
    """Harmony integration by batch_key.

    Also exports 3 UMAP plots (unintegrated, integrated 2D, integrated 3D).

    Args:
        adata (AnnData): AnnData data object. Must have 'zscore' layer and 'X_pca' obsm
            for PCA coordinates.
        batch_key (str): Column name in adata.obs containing batch identifier for integration.
        save_path (str): Directory path to save plots.
        logger (logging.Logger, optional): Optional logger object. Defaults to None.
        n_neighbors (int, optional): Number of neighbors for graph construction.
            Defaults to 20.
        n_pcs (int, optional): Number of principal components to use. Defaults to 50.

    Returns:
        AnnData: Integrated data with new embeddings:
            - obsm['X_pca_harmony']: Harmony-corrected PCA coordinates
            - obsm['X_umap_harmony_n_n']: UMAP from harmony coordinates
            - obsp['neighbors_harmony_n_n']: Neighbor graph from harmony coordinates
    """
    sc.settings.figdir = save_path
    adata.X = adata.layers["zscore"]

    if logger:
        logger.info("Integration: Run Harmony")
    sc.external.pp.harmony_integrate(adata, key=batch_key, basis="X_pca")

    if logger:
        logger.info("Integration: Compute neighbors and UMAP")

    neighbors_key = f"neighbors_harmony_{n_neighbors}_{n_pcs}"
    umap_key = f"X_umap_harmony_{n_neighbors}_{n_pcs}"

    sc.pp.neighbors(
        adata,
        use_rep="X_pca_harmony",
        key_added=neighbors_key,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
    )

    sc.tl.umap(adata, neighbors_key=neighbors_key, key_added=umap_key, n_components=2)

    _plot_integration_comparison(adata, save_path, umap_key)

    # 3D UMAP
    sc.tl.umap(
        adata,
        neighbors_key=neighbors_key,
        key_added=f"neighbors_harmony_{n_neighbors}_{n_pcs}_3D",
        n_components=3,
    )
    with mpl.rc_context({"figure.figsize": (8, 8)}):
        sc.pl.embedding(
            adata,
            basis=f"neighbors_harmony_{n_neighbors}_{n_pcs}_3D",
            color=["sample", "condition", "cell_type_mmc_raw_revised"],
            size=0.5,
            alpha=0.02,
            wspace=0.1,
            show=False,
            projection="3d",
        )
        plt.tight_layout()
        plt.savefig(
            join(save_path, "UMAP_integrated_harmony_3D.png"),
            dpi=200,
            bbox_inches="tight",
        )
        plt.close()

    return adata


def _plot_integration_comparison(adata: AnnData, save_path: str, umap_key: str) -> None:
    """Helper function to plot before/after integration comparison."""
    pt_size = 320000 / adata.shape[0]

    fig, axes = plt.subplots(
        3, 2, figsize=(20, 24), gridspec_kw={"hspace": 0.05, "wspace": 0}
    )

    fig.text(0.32, 0.9, "Unintegrated", fontsize=24, fontweight="normal", ha="center")
    fig.text(
        0.72, 0.9, "Integrated (Harmony)", fontsize=24, fontweight="normal", ha="center"
    )

    plot_configs = [
        ("sample", "Sample"),
        ("condition", "Condition"),
        ("cell_type_mmc_raw_revised", "Cell Type"),
    ]

    for i, (color_key, label) in enumerate(plot_configs):
        sc.pl.embedding(
            adata,
            basis=umap_key.replace("_harmony", ""),
            color=color_key,
            show=False,
            ax=axes[i, 0],
            size=pt_size,
            legend_loc=None,
            title="",
        )
        axes[i, 0].set_xlabel("")
        axes[i, 0].set_ylabel("")
        axes[i, 0].set_aspect("equal")

        sc.pl.embedding(
            adata,
            basis=umap_key,
            color=color_key,
            show=False,
            ax=axes[i, 1],
            size=pt_size,
            legend_loc="right margin",
            title="",
        )
        axes[i, 1].set_xlabel("")
        axes[i, 1].set_ylabel("")
        axes[i, 1].set_aspect("equal")

        axes[i, 0].text(
            -0.1,
            0.5,
            label,
            transform=axes[i, 0].transAxes,
            fontsize=20,
            fontweight="normal",
            rotation=90,
            verticalalignment="center",
            horizontalalignment="center",
        )

    plt.tight_layout()
    plt.savefig(
        join(save_path, "UMAP_integrated_harmony.png"), dpi=200, bbox_inches="tight"
    )
    plt.close()
