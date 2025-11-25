import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

from cellseg_benchmark._constants import (
    cell_type_colors,
    contamination_markers,
    selected_EC_subtypes,
)
from cellseg_benchmark.adata_utils import (
    clean_pca_umap,
    pca_umap_single,
    integration_harmony,
    normalize_counts,
)
from cellseg_benchmark.cell_annotation_utils import (
    annotate_cells_by_score,
    flag_contamination,
    score_cell_types,
)


def main():
    parser = argparse.ArgumentParser(
        description="Identify vascular subtypes."
    )
    parser.add_argument("cohort")
    parser.add_argument("seg_method")
    parser.add_argument(
        "--cell-type-col",
        default="cell_type_revised",
    )
    parser.add_argument(
        "--condition-col",
        default="condition",
    )
    parser.add_argument(
        "--batch-key",
        default="slide",
    )
    parser.add_argument(
        "--input-filename",
        default="adata_integrated.h5ad.gz",
    )
    parser.add_argument(
        "--base-path",
        type=Path,
        default=Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"),
    )
    args = parser.parse_args()

    # Logger Setup
    logger = logging.getLogger("vascular_subtyping")
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
        logger.addHandler(handler)
    
    logger.info("Starting vascular subtyping for %s / %s", args.cohort, args.seg_method)

    # Paths
    method_path = args.base_path / "analysis" / args.cohort / args.seg_method
    adatas_path = method_path / "adatas"
    plot_path = method_path / "plots" / "subtyping-vascular"

    for p in [plot_path, plot_path / "ecs-only", plot_path / "spatial", adatas_path]:
        p.mkdir(parents=True, exist_ok=True)

    logger.info("Load adata...")
    adata = sc.read_h5ad(adatas_path / args.input_filename)
    cell_type_col = args.cell_type_col

    if cell_type_col not in adata.obs.columns:
        raise KeyError(f"'{cell_type_col}' not in adata.obs")

    adata.obs[cell_type_col] = pd.Categorical(
        adata.obs[cell_type_col], categories=list(cell_type_colors.keys())
    ).remove_unused_categories()
    adata.uns[f"{cell_type_col}_colors"] = [
        cell_type_colors[c] for c in adata.obs[cell_type_col].cat.categories
    ]

    logger.info("Subset vascular cells...")
    vascular_types = ["ECs", "SMCs", "Pericytes", "VLMCs", "BAMs"]
    sub = adata[adata.obs[cell_type_col].isin(vascular_types)].copy()
    del adata

    logger.info("Filter out contaminated cells...")
    sub = flag_contamination(
        sub, contamination_markers, layer="volume_log1p_norm",
        absolute_min=1, z_threshold=2, logger=logger
    )

    with rc_context({"figure.figsize": (6, 6)}):
        sc.pl.embedding(
            sub,
            basis="X_umap_harmony_20_50",
            color=[
                "contaminated",
                "contaminated_Neuronal",
                "contaminated_Oligodendrocytes",
                "contaminated_Astrocytes",
                "contaminated_Immune",
            ],
            size=50000 / sub.shape[0],
            show=False,
        )
        plt.savefig(plot_path / "UMAP_cell_type_contaminations.png", dpi=200)
        plt.close()

    sub = sub[~sub.obs["contaminated"]].copy()

    with rc_context({"figure.figsize": (6,6)}):
        sc.pl.embedding(
            sub, 
            basis="X_umap_harmony_20_50", 
            color=cell_type_col, 
            show=False,
        )
        plt.savefig(
            plot_path / "UMAP_integrated_subset.png", dpi=200, bbox_inches="tight"
        )
        plt.close()
    
    logger.info("Re-process vascular subset...")
    sub = clean_pca_umap(sub, logger=logger)
    sub.X = sub.layers["counts"]
    assert np.issubdtype(sub.X.dtype, np.integer)
    sub = normalize_counts( # do not trim cells again in subset
        sub, save_path=None, seg_method=args.seg_method, trim_outliers=False, logger=logger
    )
    sub = pca_umap_single(
        sub, n_neighbors=10, n_pcs=20, save_path=None, logger=logger
    )
    sub = integration_harmony(
        sub,
        batch_key=args.batch_key,
        n_neighbors=10,
        n_pcs=20,
        save_path=None,
        logger=logger,
    )

    logger.info("Score EC zonation subtypes...")
    # Step 1: score marker genes
    sub = score_cell_types(
        sub, selected_EC_subtypes, top_n_genes=10, layer="volume_log1p_norm"
    )
    # Step 2: assign each cell by max scoring subtype
    sub = annotate_cells_by_score(
        sub, selected_EC_subtypes, out_col="ec_zonation", score_threshold=0.15
    )

    sub.obs["ec_zonation"] = pd.Categorical(
        sub.obs["ec_zonation"], categories=list(cell_type_colors.keys())
    ).remove_unused_categories()
    sub.uns["ec_zonation_colors"] = [
        cell_type_colors[ct] for ct in sub.obs["ec_zonation"].cat.categories
    ]

    logger.info("Re-process EC-only subset...")
    ec = sub[sub.obs["ec_zonation"].isin(["aECs", "capECs", "vECs"])].copy()
    logger.info(f"Cells removed (non-ECs): {sub.shape[0]-ec.shape[0]}")
    logger.info(f"Remaining ECs: {ec.shape[0]}")
    
    ec = clean_pca_umap(ec, logger=logger)
    ec.X = ec.layers["counts"]
    assert np.issubdtype(ec.X.dtype, np.integer)
    ec = normalize_counts( # do not trim cells again in subset
        ec, save_path=None, seg_method=args.seg_method, trim_outliers=False, logger=logger
    )
    ec = pca_umap_single(
        ec, n_neighbors=10, n_pcs=20, save_path=None, logger=logger
    )
    ec = integration_harmony(
        ec,
        batch_key=args.batch_key,
        n_neighbors=10,
        n_pcs=20,
        save_path=None,
        logger=logger,
    )

    with rc_context({"figure.figsize": (6,6)}):
        sc.pl.embedding(
            ec,
            basis="X_umap_harmony_10_20",
            color=[cell_type_col, "ec_zonation", "score_aECs", "score_capECs", "score_vECs"], 
            cmap="inferno",
            ncols=5,
            wspace=0.25,
            size=320000 / ec.shape[0],
            show=False,
        )
        plt.savefig(
            plot_path / "ecs-only" / "UMAP_ECs.png", dpi=150, bbox_inches="tight"
        )
        plt.close()
        
    # run 3D UMAP (computation affects 2D components, therefore save separately)
    sc.pp.neighbors(ec, n_neighbors=10, n_pcs=20)
    sc.tl.umap(ec, neighbors_key="neighbors_harmony_10_20", key_added="X_umap_harmony_10_20_3D", n_components=3)
    with rc_context({'figure.figsize': (6,6)}):
        sc.pl.embedding(ec, 
                        basis="X_umap_harmony_10_20_3D", 
                        color=[cell_type_col, "sample", "ec_zonation"], 
                        alpha=0.2, 
                        wspace=0.1, 
                        projection='3d', 
                        show=False,
                       )
        plt.savefig(
            plot_path / "ecs-only" / "UMAP_ECs_3D.png", dpi=150, bbox_inches="tight"
        )
        plt.close()
        
    with rc_context({'figure.figsize': (6,6)}):
        sc.pl.embedding(ec, 
                        basis="X_umap_harmony_10_20", 
                        color=['n_counts','n_genes', 'volume_final'], 
                        cmap='turbo',
                        size=320000 / ec.shape[0], 
                        show=False,
                       )
        plt.savefig(
            plot_path / "ecs-only" / "UMAP_ECs_metadata.png", dpi=100, bbox_inches="tight"
        )
        plt.close()

    sc.external.tl.trimap(ec)
    with rc_context({"figure.figsize": (6,6)}):
        sc.external.pl.trimap(ec,
                              color=[cell_type_col, "ec_zonation", "score_aECs", "score_capECs", "score_vECs"], 
                              cmap="inferno",
                              ncols=5,
                              wspace=0.25
                              size=320000 / ec.shape[0], 
                              show=False, 
                             )
        plt.savefig(
            plot_path / "ecs-only" / "TRIMAP_ECs.png", dpi=150, bbox_inches="tight"
        )
        plt.close()
    
    sc.external.tl.phate(ec)
    with rc_context({"figure.figsize": (6,6)}):
        sc.external.pl.phate(ec, ncols=5,
                             color=[cell_type_col, "ec_zonation", "score_aECs", "score_capECs", "score_vECs"], 
                             cmap="inferno",
                             ncols=5,
                             wspace=0.25
                             size=320000 / ec.shape[0], 
                             show=False, 
                            )
        plt.savefig(
            plot_path / "ecs-only" / "PHATE_ECs.png", dpi=150, bbox_inches="tight"
        )
        plt.close()

    logger.info("Summaries by condition...")
    cond = args.condition_col
    if cond in ec.obs.columns:
        df_counts = pd.crosstab(ec.obs[cond], ec.obs["ec_zonation"], margins=True)
        df_frac = pd.crosstab(ec.obs[cond], ec.obs["ec_zonation"], normalize="index", margins=True).mul(100).round(1)
        df_frac_sample = pd.crosstab(
            [ec.obs[cond], ec.obs["sample"]],
            ec.obs["ec_zonation"],
            normalize='index',
            margins=True
        ).mul(100).round(1)
        logger.info("EC subtype counts:\n%s", df_counts)
        logger.info("EC subtype fractions:\n%s", df_frac)
        logger.info("EC subtype fractions (per sample):\n%s", df_frac_sample)

    logger.info("Merge EC zonation subtypes back into vascular set...")
    merged = sub.copy()
    if "ec_zonation" in merged.obs.columns:
        merged.obs = merged.obs.drop(columns="ec_zonation")
    
    merged.obs = merged.obs.merge(
        ec.obs[["ec_zonation"]], left_index=True, right_index=True, how="left"
    )
    merged.obs["cell_type_incl_zonation"] = (
        merged.obs["ec_zonation"]
        .astype("object")
        .fillna(merged.obs[cell_type_col].astype("object"))
    )
    
    # set "otherECs" for specific cell types without ec annotation. these were previously removed as contamination.
    missed_ECs_mask = merged.obs["ec_zonation"].isna() & merged.obs[
        cell_type_col
    ].isin(["ECs"])
    merged.obs.loc[missed_ECs_mask, "cell_type_incl_zonation"] = "otherECs"
    
    merged.obs["cell_type_incl_zonation"] = merged.obs["cell_type_incl_zonation"].astype("category")
    merged.obs.drop(columns="ec_zonation", inplace=True)
    
    # Set ordered categories and color map
    merged.obs["cell_type_incl_zonation"] = pd.Categorical(
        merged.obs["cell_type_incl_zonation"], categories=list(cell_type_colors.keys())
    ).remove_unused_categories()
    
    merged.uns["cell_type_incl_zonation_colors"] = [
        cell_type_colors[ct] for ct in merged.obs["cell_type_incl_zonation"].cat.categories
    ]

    with rc_context({"figure.figsize": (6,6)}):
        sc.pl.embedding(
            merged,
            basis="X_umap_harmony_10_20",
            color=[cell_type_col, "cell_type_incl_zonation", "score_aECs", "score_capECs", "score_vECs"], 
            cmap="inferno",
            ncols=5,
            wspace=0.25
            size=320000 / ec.shape[0], 
            show=False, 
           )
        plt.savefig(
            plot_path / "UMAP_integrated_EC_scoring.png", dpi=150, bbox_inches="tight", 
        )
        plt.close()
    
    logger.info(f"Exporting spatial plots to: {plot_path/"spatial"}...")
    if "sample" in merged.obs.columns:
        for s in merged.obs["sample"].unique():
            with rc_context({"figure.figsize": (10, 10)}):
                sc.pl.embedding(
                    merged[merged.obs["sample"] == s],
                    basis="spatial_microns",
                    color="cell_type_incl_zonation",
                    size=9,
                    show=False,
                )
                plt.savefig(plot_path / "spatial" / f"Spatial_{s}.png", dpi=500)
                plt.close()
    
    out = adatas_path / "adata_vascular_subset.h5ad.gz"
    logger.info("Saving to %s", out)
    if "fov" in merged.obs.columns:
        merged.obs["fov"] = merged.obs["fov"].astype(str)

    merged.write(out, compression="gzip")
    logger.info("Done.")


if __name__ == "__main__":
    main()