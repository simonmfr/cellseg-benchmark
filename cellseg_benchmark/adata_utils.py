import logging

from matplotlib import pyplot as plt
from spatialdata import SpatialData
from anndata import AnnData, concat
from typing import List, Tuple
from tqdm import tqdm
import json
import seaborn as sns
from os.path import join
import scanpy as sc
import numpy as np
import squidpy as sq
import matplotlib as mpl
import scipy.sparse as sp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from cellseg_benchmark.cell_annotation_utils import cell_type_colors

def merge_adatas(sdatas: List[Tuple[str, SpatialData]], key: str, logger: logging.Logger=None, do_qc=False, save_path=None) -> AnnData:
    adatas = []
    with open(
            "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
    ) as f:
        paths = json.load(f)

    y_limits = [0,0,0,0]
    if logger:
        logger.info(f"Merging adatas of {key}")
    for name, sdata in tqdm(sdatas):
        if f"adata_{key}" not in sdata.tables.keys():
            if logger:
                logger.warning(f"Skipping {name}. No such key: {key}")
            continue
        adata = sdata[f"adata_{key}"]
        samples = name.split("_")
        adata.obs['cohort'] = samples[0]
        adata.obs['slide'] = samples[1]
        adata.obs['region'] = samples[2]
        adata.obs['condition'] = paths[name].split("/")[-2].split("-")[-2]
        if isinstance(adata.X, np.ndarray):
            adata.X = sp.csr_matrix(adata.X, dtype=np.float32)
        adata.obs['n_counts'] = adata.X.sum(axis=1)
        adata.obs['n_genes'] = adata.X.count_nonzero(axis=1)
        adata.obs['sample'] = name
        adatas.append(adata)

        if do_qc:
            y_limits[0] = max(y_limits[0], max(np.histogram(adata.obs['n_counts'], bins=60)[0]))
            y_limits[1] = max(y_limits[1], max(np.histogram(adata.obs['n_counts'][adata.obs["n_counts"] < 1000], bins=60)[0]))
            y_limits[2] = max(y_limits[2], max(np.histogram(adata.obs["n_genes"], bins=60)[0]))
            y_limits[3] = max(y_limits[3], max(np.histogram(adata.obs["volume"], bins=100)[0]))

    adata = concat(adatas, join="outer", merge="first")
    for i in set("cell_id", "cell_id_x", "cell_id_y") & set(adata.obs.columns):
        del adata.obs[i]
    adata.obs["cell_type_mmc_raw_revised"] = pd.Categorical(
        adata.obs["cell_type_mmc_raw_revised"], categories=list(cell_type_colors.keys())
    )
    adata.obs["cell_type_mmc_raw_revised"] = adata.obs["cell_type_mmc_raw_revised"].cat.remove_unused_categories()

    # set colors in adata.uns
    colors = [cell_type_colors[cat] for cat in adata.obs["cell_type_mmc_raw_revised"].cat.categories]
    adata.uns["cell_type_mmc_raw_revised_colors"] = colors
    adata.obs_names_make_unique()
    if logger:
        logger.info(f"{key}: # of cells: {len(adata)}, # of samples: {len(adata.obs.sample.unique())}")
    if do_qc:
        plot_qc(adata, save_path, y_limits, logger)
    return adata

def merge_adatas_deb(sdatas: List[Tuple[str, AnnData]], do_qc=False, save_path=None, logger=None) -> AnnData:
    adatas = []
    paths = {
  "foxf2_s1_r0": "/dss/dssfs03/pn52re/pn52re-dss-0000/202402-Foxf2/merfish_output/20240229_Foxf2-Slide01-cp-WT-ECKO/region_0-ECKO000",
  "foxf2_s1_r1": "/dss/dssfs03/pn52re/pn52re-dss-0000/202402-Foxf2/merfish_output/20240229_Foxf2-Slide01-cp-WT-ECKO/region_1-WT000"
    }

    y_limits = [0,0,0,0]

    for name, adata in tqdm(sdatas):
        samples = name.split("_")
        adata.obs['cohort'] = samples[0]
        adata.obs['slide'] = samples[1]
        adata.obs['region'] = samples[2]
        adata.obs['condition'] = paths[name].split("/")[-2].split("-")[-2]
        adata.obs['n_counts'] = adata.X.sum(axis=1)
        adata.obs['n_genes'] = adata.X.count_nonzero(axis=1)
        adata.obs['sample'] = name
        adata.obs['volume'] = 7*adata.obs.area/(10*10)
        adatas.append(adata)

        if do_qc:
            y_limits[0] = max(y_limits[0], max(np.histogram(adata.obs['n_counts'], bins=60)[0]))
            y_limits[1] = max(y_limits[1], max(np.histogram(adata.obs['n_counts'][adata.obs["n_counts"] < 1000], bins=60)[0]))
            y_limits[2] = max(y_limits[2], max(np.histogram(adata.obs["n_genes"], bins=60)[0]))
            y_limits[3] = max(y_limits[3], max(np.histogram(adata.obs["volume"], bins=100)[0]))

    adata = concat(adatas, join="outer", merge="first")
    adata.obs_names_make_unique()
    if do_qc:
        plot_qc(adata, save_dir=save_path, y_limits=y_limits, logger=logger)
    return adata

def plot_qc(adata: AnnData, save_dir, y_limits, logger) -> None:
    if logger:
        logger.info(f"Plotting QC results")
    sample_names = adata.obs.sample.unique()

    #=====================================General Stats=================================================================
    fig, axs = plt.subplots(len(sample_names), 4, figsize=(16, len(sample_names) * 4), gridspec_kw={'wspace': 0.4}, squeeze=False)
    for i, name in enumerate(sample_names):
        adata_tmp = adata[adata.obs["sample"] == name]

        for ax, data, bins, y_limit in zip(
                axs[i,:],
                [adata_tmp.obs["n_counts"],
                 adata_tmp.obs["n_counts"][adata_tmp.obs["n_counts"] < 1000],
                 adata_tmp.obs["n_genes"],
                 adata_tmp.obs["volume"]],
                [60, 60, 60, 100],
                y_limits
        ):
            sns.histplot(data, kde=False, bins=bins, ax=ax, zorder=2, alpha=1)
            ax.grid(True, zorder=1)
            ax.set_ylim(0, y_limit)  # Set consistent y-axis limit for each column
            for patch in ax.patches:
                patch.set_edgecolor('none')  # Remove border lines

    pad = 5
    cols = ["Number of Counts per Cell", "Number of Counts per Cell (until n=1000)", "Number of Genes per Cell", "Volume"]
    rows = sample_names

    for ax, col in zip(axs[0,:], cols):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                    xycoords='axes fraction', textcoords='offset points',
                    size='large', ha='center', va='baseline')

    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center', rotation='vertical')

    fig.tight_layout()
    fig.savefig(join(save_dir, "qc_general_stats_cells.png"))
    plt.close()

    #===================================Cell Qualities==================================================================
    fig, axs = plt.subplots(len(sample_names), 1, figsize=(9, len(sample_names) * 8),
                            gridspec_kw={'wspace': 0.4})
    for ax, name in zip(axs, sample_names):
        adata_tmp = adata[adata.obs["sample"] == name]
        sc.pl.scatter(adata_tmp, x='n_counts', y='n_genes', color="volume", ax=ax, show=False, legend_loc='none')

        default_cax = fig.axes[-1]
        if default_cax is not ax:
            fig.delaxes(default_cax)

        im = ax.collections[-1]
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.1)
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label('volume')

    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center', rotation='vertical')

    fig.tight_layout()
    fig.savefig(join(save_dir, "qc_cell_qualities.png"))
    plt.close()

    #===============================General Stats Genes=================================================================
    fig, axs = plt.subplots(len(sample_names), 3, figsize=(18,len(sample_names) * 5),
                            gridspec_kw={'wspace': 0.4})
    for i, name in enumerate(sample_names):
        adata_tmp = adata[adata.obs["sample"] == name]
        adata_tmp.var["n_counts"] = adata_tmp.X.sum(axis=0).T
        log10_counts = np.log10(adata_tmp.var["n_counts"] + 1)

        sns.violinplot(y=log10_counts, ax=axs[i, 0], zorder=2)
        axs[i, 0].set_ylabel("log10(n_counts+1)")
        axs[i, 0].grid(True, alpha=0.5, zorder=0)

        sns.histplot(log10_counts, kde=True, bins=25, ax=axs[i,1], zorder=2)
        axs[i,1].set_xlabel("log10(n_counts+1)")
        axs[i,1].grid(True, alpha=0.5, zorder=0)

        sns.ecdfplot(adata_tmp.var["n_counts"], log_scale=True, ax=axs[i,2], zorder=2)
        axs[i,2].grid(True, alpha=0.5, zorder=0)

    for ax, row in zip(axs[:,0], rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center', rotation='vertical')

    fig.tight_layout()
    fig.savefig(join(save_dir, "qc_general_stats_genes.png"))
    plt.close()

    #===============================Highly expressed genes==============================================================
    fig, axs = plt.subplots(len(sample_names), 1, figsize=(12, len(sample_names) * 4),
                            gridspec_kw={'wspace': 0.4})
    for ax, name in zip(axs, sample_names):
        adata_tmp = adata[adata.obs["sample"] == name]
        sc.pl.highest_expr_genes(adata_tmp, n_top=20, ax=ax, show=False)

    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center', rotation='vertical')

    fig.tight_layout()
    fig.savefig(join(save_dir, "qc_highes_expressed_genes.png"))
    plt.close()

def filter_cells(adata, save_path, min_counts=25, min_genes=5, logger=None):
    sample_names = adata.obs["sample"].unique()

    fig, axs = plt.subplots(len(sample_names), 1, figsize=(12, len(sample_names)*9), gridspec_kw={'wspace': 0.4})
    for ax, name in zip(axs, sample_names):
        adata_tmp = adata[adata.obs["sample"] == name]
        adata_tmp.obs["cell_outlier"] = (adata_tmp.obs["n_counts"] < min_counts) | (
        adata_tmp.obs["n_genes"] < min_genes)
        adata_tmp.obs["cell_outlier"] = adata_tmp.obs["cell_outlier"].astype('category')
        sq.pl.spatial_scatter(
            adata_tmp, ax=ax, shape=None, color="cell_outlier", size=0.125, library_id="spatial", figsize=(8, 8),
            palette=mpl.colors.ListedColormap(["lightgrey", "red"])
        )

    pad = 5
    rows = sample_names
    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center', rotation='vertical')
    fig.savefig(join(save_path, f"filter_outlier.png"))
    plt.close()

    if logger:
        logger.info(f"# cells before filtering: {adata.n_obs}")
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if logger:
        logger.info(f"# cells after filtering: {adata.n_obs}")
        n = (adata.obs["volume"] < 100).sum()
        logger.info(f"{n} cells with volume <100 um3")
    metric = "volume"
    n = 3
    adata.obs["outlier"] = (
            (adata.obs[metric] > n * np.median(adata.obs[metric])) |
            (adata.obs[metric] < 100)
    ).astype('category')
    if logger:
        logger.info(adata.obs.outlier.value_counts())
    fig, axs = plt.subplots(len(sample_names), 1, figsize=(12, len(sample_names)*9), gridspec_kw={'wspace': 0.4})
    for ax,name in zip(axs, sample_names):
        adata_tmp = adata[adata.obs["sample"] == name]
        sq.pl.spatial_scatter(
            adata_tmp, ax=ax, shape=None, color="outlier", size=0.125, library_id="spatial", figsize=(8, 8),
            palette=mpl.colors.ListedColormap(["lightgrey", "red"]), save=join(save_path, f"volume_outlier_{name}.png")
        )

    pad = 5
    rows = sample_names
    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center', rotation='vertical')
    fig.savefig(join(save_path, f"filter_volume_outlier.png"))
    plt.close()
    adata = adata[adata.obs['outlier'] == False]
    return adata

def filter_genes(adata, save_path, logger=None):
    if logger:
        logger.info(f"# genes before filtering: {adata.n_vars}")
    sc.pp.filter_genes(adata, min_cells=10)
    if logger:
        logger.info(f"# genes after filtering: {adata.n_vars}")

    fig, axs = plt.subplots(1, 4, figsize=(18, 4), gridspec_kw={'wspace': 0.4})
    for ax, data, bins in zip(
            axs,
            [adata.obs["n_counts"],
             adata.obs["n_counts"][adata.obs["n_counts"] < 1000],
             adata.obs["n_genes"],
             adata.obs["volume"]],
            [100, 100, 100, 100]  # bins
    ):
        sns.histplot(data, kde=False, bins=bins, ax=ax, zorder=2, alpha=1)
        ax.grid(True, zorder=1)
        for patch in ax.patches:
            patch.set_edgecolor('none')
    fig.savefig(join(save_path, "filter_general_stats_merged.png"))
    plt.close()
    return adata

def normalize(adata, save_path, logger=None):
    def normalize_after_volume(adata, save_path, logger=None):
        if logger:
            logger.info(f"normalize after volume")
        adata.layers["counts"] = adata.X
        volumes = adata.obs['volume'].values
        sparse_matrix = sp.csr_matrix(adata.X, dtype=np.float32)
        normalized_matrix = sparse_matrix.multiply(1 / volumes[:, None])
        normalized_matrix2 = sp.csr_matrix(normalized_matrix, dtype=np.float32)
        adata.layers["volume_norm"] = normalized_matrix2
        del normalized_matrix, normalized_matrix2
        return adata

    def remove_outliers(adata, save_path, logger=None):
        if logger:
            logger.info(f"remove outliers")
        row_sums = adata.layers["volume_norm"].sum(axis=1).A1
        p1 = np.percentile(row_sums, 1)
        p2 = np.percentile(row_sums, 2)
        p98 = np.percentile(row_sums, 98)
        p99 = np.percentile(row_sums, 99)
        sns.histplot(row_sums)
        plt.axvline(p1, color='red', lw=1)
        plt.axvline(p99, color='red', lw=1)  # less strict
        plt.axvline(p2, color='orange', lw=1)
        plt.axvline(p98, color='orange', lw=1)  # used by Allen et al.
        plt.savefig(join(save_path, "outlier_counts.png"))
        plt.close()
        mask = (row_sums > p1) & (row_sums < p99)
        adata = adata[mask]
        if logger:
            logger.info(f"# cells after filtering: {adata.n_obs}")
        return adata

    def z_score(adata, save_path, logger=None):
        if logger:
            logger.info(f"calculate z-score")
        adata.layers["zscore"] = sc.pp.scale(adata.layers["volume_log1p_norm"], zero_center=True, max_value=None)
        if logger:
            logger.info(f"normalize after lib size")
        scales_counts = sc.pp.normalize_total(adata, layer="counts", target_sum=None, inplace=False)
        adata.layers["librarysize_log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

        fig, axs = plt.subplots(1, 5, figsize=(25, 4), gridspec_kw={'wspace': 0.4})
        titles = ['Counts', 'Normalized (by volume)', 'Log-normalized (by volume)', "Z-score (after volume log-norm)", 'Log-normalized (by libsize)']
        for ax, data, bins, title in zip(
                axs,
                [adata.layers["counts"].sum(1),
                 adata.layers["volume_norm"].sum(1),
                 adata.layers["volume_log1p_norm"].sum(1),
                 adata.layers["zscore"].sum(1),
                 adata.layers["librarysize_log1p_norm"].sum(1)],
                [100, 10, 100, 100, 100],
                titles
        ):
            sns.histplot(data, kde=False, bins=bins, ax=ax, zorder=2, alpha=1)
            ax.set_title(title)
            ax.grid(True)
            legend = ax.get_legend()
            if legend is not None:
                legend.remove()
            for patch in ax.patches:
                patch.set_edgecolor('none')
        plt.savefig(join(save_path, "normalize_z_score.png"))
        plt.close()
        return adata

    def scaling(adata, scaling_factor, logger=None):
        if logger:
            logger.info(f"scaling to {scaling_factor}")
        normalized_matrix = adata.layers["volume_norm"]

        # Scale so that the sum of values per cell (per row)  is equal to 250, as done by Allen et al. 2024 Cell
        row_sums = normalized_matrix.sum(axis=1).A1
        scaling_factors = scaling_factor / row_sums
        normalized_matrix2 = normalized_matrix.multiply(scaling_factors[:, None])

        normalized_matrix2 = sp.csr_matrix(normalized_matrix2, dtype=np.float32)
        adata.layers["volume_norm"] = normalized_matrix2
        assert (round(normalized_matrix2.sum(axis=1).A1.mean()) == 250)
        del normalized_matrix2, normalized_matrix
        return adata

    adata = normalize_after_volume(adata, save_path, logger=logger)
    adata = remove_outliers(adata, save_path, logger=logger)
    adata = scaling(adata, scaling_factor=250, logger=logger)
    adata.layers[("volume_log1p_norm")] = sc.pp.log1p(adata.layers["volume_norm"], copy=True)
    adata = z_score(adata, save_path, logger=logger)
    return adata

def dimensionality_reduction(adata, save_path, logger=None):
    adata.X = adata.layers["zscore"]
    sc.settings.figdir = save_path
    if logger:
        logger.info("Dimensionality reduction: PCA")
    sc.pp.pca(adata, svd_solver="arpack")
    fig, axs = plt.subplots(1, 2, figsize=(20, 9), gridspec_kw={'wspace': 0.4})
    sc.pl.pca_scatter(adata, color="n_counts", ax=axs[0])
    var_ratio = adata.uns['pca']['variance_ratio']

    # plot on your second subplot
    axs[1].plot(
        np.arange(1, len(var_ratio) + 1),
        var_ratio,
        marker='o',
        linestyle='-'
    )
    axs[1].set_yscale('log')
    axs[1].set_xlabel('Principal component')
    axs[1].set_ylabel('Variance ratio')
    axs[1].set_title('PCA variance ratio')
    fig.savefig(join(save_path, "PCA.png"))
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(30, 8), gridspec_kw={'wspace': 0.4})
    if logger:
        logger.info("Dimensionality reduction: neighbors and umap with n=20 PCA")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    with plt.style.context('default'):
        with mpl.rc_context({'figure.figsize': (8, 8)}):
            sc.pl.embedding(adata, ax=axs[0], basis="X_umap", color="n_counts", title="UMAP: n = 20", show=False)

    if logger:
        logger.info("Dimensionality reduction: neighbors and umap with n=40 PCA")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    with plt.style.context('default'):
        with mpl.rc_context({'figure.figsize': (8, 8)}):
            sc.pl.embedding(adata, ax=axs[1], basis="X_umap", color="n_counts", title="UMAP: n = 40", show=False)

    if logger:
        logger.info("Dimensionality reduction: neighbors and umap with n=50 PCA")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
    sc.tl.umap(adata)
    with plt.style.context('default'):
        with mpl.rc_context({'figure.figsize': (8, 8)}):
            sc.pl.embedding(adata, ax=axs[2], basis="X_umap", color="n_counts", title="UMAP: n = 50", show=False)

    plt.savefig(join(save_path, "UMAP_unintegrated.png"))
    plt.close()
    return adata

def integration_harmony(adata, batch_key, save_path, logger=None):
    sc.settings.figdir = save_path
    adata.X = adata.layers["zscore"]
    if logger:
        logger.info("Integration harmony: PCA")
    sc.pp.pca(adata, svd_solver="arpack", key_added="X_pca_new")
    if logger:
        logger.info("Integration harmony: harmony")
    sc.external.pp.harmony_integrate(adata, key=batch_key, basis="X_pca_new")
    del adata.obsm["X_pca_new"]

    if logger:
        logger.info("Integration harmony: neighbors")
    sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=50, key_added="neighbors_harmony")
    if logger:
        logger.info("Integration harmony: umap")
    sc.tl.umap(adata, neighbors_key="neighbors_harmony", key_added="X_umap_harmony")

    pt_size_umap = 220000 / adata.shape[0]
    fig, axs = plt.subplots(3, 2, figsize=(30, 35), gridspec_kw={'wspace': 0.4})
    sc.pl.embedding(adata, basis="X_umap", color='sample', show=False, ax=axs[0,0], title="No Integration", size=pt_size_umap, legend_loc=None)
    axs[0,0].set_aspect('equal')
    sc.pl.embedding(adata, basis="X_umap_harmony", color='sample', show=False, ax=axs[0,1], title="Harmony Integration", size=pt_size_umap)
    axs[0,1].set_aspect('equal')
    sc.pl.embedding(adata, basis="X_umap", color='slide', show=False, ax=axs[1,0], title="No Integration", size=pt_size_umap, legend_loc=None)
    axs[1,0].set_aspect('equal')
    sc.pl.embedding(adata, basis="X_umap_harmony", color='slide', show=False, ax=axs[1,1], title="Harmony Integration", size=pt_size_umap)
    axs[1,1].set_aspect('equal')
    sc.pl.embedding(adata, basis="X_umap", color='cell_type_mmc_raw_revised', show=False, ax=axs[2,0], title="No Integration", size=pt_size_umap, legend_loc=None)
    axs[2,0].set_aspect('equal')
    sc.pl.embedding(adata, basis="X_umap_harmony", color='cell_type_mmc_raw_revised', show=False, ax=axs[2,1], title="Harmony Integration", size=pt_size_umap)
    axs[2,1].set_aspect('equal')

    pad = 5
    rows = ["Sample", "Slide", "Cell Type"]

    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center', rotation='vertical')
    plt.savefig(join(save_path, "UMAP_integration.png"))
    plt.close()
    return adata
