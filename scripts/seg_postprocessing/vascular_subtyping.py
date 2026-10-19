#!/usr/bin/env python
import argparse
import logging
import pathlib
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))

import cellseg_benchmark._constants as _constants
import cellseg_benchmark._markers as _markers
import cellseg_benchmark.adata_utils as adata_utils
import cellseg_benchmark.cell_annotation_utils as cell_annotation_utils

logger = logging.getLogger("vascular_subtyping")
logger.setLevel(logging.INFO)
_handler = logging.StreamHandler()
_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(_handler)


def plot_embedding(adata, basis, color, path, ncols=5):
    with plt.rc_context({"figure.figsize": (6, 6)}):
        sc.pl.embedding(
            adata,
            basis=basis,
            color=color,
            cmap="inferno",
            ncols=ncols,
            wspace=0.25,
            size=320000 / adata.shape[0],
            show=False,
        )
        plt.savefig(path, dpi=150, bbox_inches="tight")
        plt.close()


def main():
    parser = argparse.ArgumentParser(description="Identify vascular subtypes.")
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
        type=pathlib.Path,
        default=pathlib.Path(_constants.BASE_PATH),
    )
    args = parser.parse_args()

    logger.info("Starting vascular subtyping for %s / %s", args.cohort, args.seg_method)

    method_path = args.base_path / "analysis" / args.cohort / args.seg_method
    adatas_path = method_path / "adatas"
    plot_path = method_path / "plots" / "subtyping-vascular"
    out = adatas_path / "adata_vascular_subset.h5ad.gz"

    out.unlink(missing_ok=True)
    shutil.rmtree(plot_path, ignore_errors=True)
    (plot_path / "spatial").mkdir(parents=True)

    logger.info("Load adata...")
    adata = sc.read_h5ad(adatas_path / args.input_filename)
    cell_type_col = args.cell_type_col

    if cell_type_col not in adata.obs.columns:
        raise KeyError(f"'{cell_type_col}' not in adata.obs")

    missing = [
        k
        for k, g in _markers.selected_EC_subtypes.items()
        if not set(g) & set(adata.var_names)
    ]
    if missing:
        raise ValueError(f"No zonation markers on panel for {missing}")

    adata.obs[cell_type_col] = pd.Categorical(
        adata.obs[cell_type_col], categories=list(_constants.cell_type_colors.keys())
    ).remove_unused_categories()
    adata.uns[f"{cell_type_col}_colors"] = [
        _constants.cell_type_colors[c] for c in adata.obs[cell_type_col].cat.categories
    ]

    logger.info("Subset vascular cells...")
    vascular_types = ["ECs", "SMCs", "Pericytes", "VLMCs", "BAMs"]
    sub = adata[adata.obs[cell_type_col].isin(vascular_types)].copy()
    del adata

    logger.info("Filter out contaminated cells...")
    sub = cell_annotation_utils.flag_contamination(
        sub,
        _markers.contamination_markers,
        layer="volume_log1p_norm",
        absolute_min=1,
        z_threshold=2,
        logger=logger,
    )

    plot_embedding(
        sub,
        "X_umap_harmony_20_50",
        [cell_type_col, "contaminated"]
        + [f"contaminated_{k}" for k in _markers.contamination_markers],
        plot_path / "vascular_contamination_umap.png",
        ncols=6,
    )
    sub = sub[~sub.obs["contaminated"]].copy()

    logger.info("Re-process vascular subset...")
    sub = adata_utils.clean_pca_umap(sub, logger=None)
    sub.X = sub.layers["counts"]
    assert np.issubdtype(sub.X.dtype, np.integer)
    sub = adata_utils.normalize_counts(  # do not trim cells again in subset
        sub,
        save_path=None,
        seg_method=args.seg_method,
        trim_outliers=False,
        logger=logger,
    )
    sub = adata_utils.pca_umap_single(
        sub, n_neighbors=10, n_pcs=20, save_path=None, logger=logger
    )
    sub = adata_utils.integration_harmony(
        sub,
        batch_key=args.batch_key,
        n_neighbors=10,
        n_pcs=20,
        save_path=None,
        logger=logger,
    )

    logger.info("Score EC zonation subtypes...")
    sub = cell_annotation_utils.score_cell_types(
        sub, _markers.selected_EC_subtypes, top_n_genes=10, layer="volume_log1p_norm"
    )
    sub = cell_annotation_utils.annotate_cells_by_score(
        sub,
        _markers.selected_EC_subtypes,
        out_col="ec_zonation",
        score_threshold=0.15,
    )

    logger.info("Gate non-EC rescue by pan-EC score...")
    top100 = pd.read_csv(
        args.base_path
        / "misc/scRNAseq_ref_ABCAtlas_Yao2023Nature/marker_genes_df/20260927_cell_type_markers_top100.csv"
    )
    pan_ec = [g for g in top100["ECs"] if g in sub.var_names][:10]
    sub = cell_annotation_utils.score_cell_types(
        sub, {"panEC": pan_ec}, top_n_genes=10, layer="volume_log1p_norm"
    )
    is_ec = sub.obs[cell_type_col] == "ECs"
    cutoff = sub.obs.loc[is_ec, "score_panEC"].quantile(0.05)
    ec_like = is_ec | (sub.obs["score_panEC"] >= cutoff)
    logger.info(f"Gated non-ECs: {(~ec_like).sum()}")
    sub.obs["ec_zonation"] = sub.obs["ec_zonation"].astype("object").where(ec_like)

    sub.obs["ec_zonation"] = pd.Categorical(
        sub.obs["ec_zonation"], categories=list(_constants.cell_type_colors.keys())
    ).remove_unused_categories()
    sub.uns["ec_zonation_colors"] = [
        _constants.cell_type_colors[ct] for ct in sub.obs["ec_zonation"].cat.categories
    ]

    logger.info("Re-process EC-only subset...")
    ec = sub[sub.obs["ec_zonation"].isin(["aECs", "capECs", "vECs"])].copy()
    logger.info(f"Cells removed (non-ECs): {sub.shape[0] - ec.shape[0]}")
    logger.info(f"Remaining ECs: {ec.shape[0]}")

    ec = adata_utils.clean_pca_umap(ec, logger=None)
    ec.X = ec.layers["counts"]
    assert np.issubdtype(ec.X.dtype, np.integer)
    ec = adata_utils.normalize_counts(  # do not trim cells again in subset
        ec,
        save_path=None,
        seg_method=args.seg_method,
        trim_outliers=False,
        logger=logger,
    )
    ec = adata_utils.pca_umap_single(
        ec, n_neighbors=10, n_pcs=20, save_path=None, logger=logger
    )
    ec = adata_utils.integration_harmony(
        ec,
        batch_key=args.batch_key,
        n_neighbors=10,
        n_pcs=20,
        save_path=None,
        logger=logger,
    )

    zonation = [
        cell_type_col,
        "ec_zonation",
        "score_aECs",
        "score_capECs",
        "score_vECs",
    ]
    plot_embedding(ec, "X_umap_harmony_10_20", zonation, plot_path / "ecs_umap.png")
    plot_embedding(
        ec,
        "X_umap_harmony_10_20",
        ["n_counts", "n_genes", "volume_final"],
        plot_path / "ecs_umap_qc.png",
    )
    logger.info("Run Trimap...")
    sc.external.tl.trimap(ec)
    plot_embedding(ec, "X_trimap", zonation, plot_path / "ecs_trimap.png")
    logger.info("Run Phate...")
    sc.external.tl.phate(ec)
    plot_embedding(ec, "X_phate", zonation, plot_path / "ecs_phate.png")

    logger.info("Summaries by condition...")
    cond = args.condition_col
    if cond in ec.obs.columns:
        df_counts = pd.crosstab(ec.obs[cond], ec.obs["ec_zonation"], margins=True)
        df_frac = (
            pd.crosstab(
                ec.obs[cond], ec.obs["ec_zonation"], normalize="index", margins=True
            )
            .mul(100)
            .round(1)
        )
        df_frac_sample = (
            pd.crosstab(
                [ec.obs[cond], ec.obs["sample"]],
                ec.obs["ec_zonation"],
                normalize="index",
                margins=True,
            )
            .mul(100)
            .round(1)
        )
        logger.info("EC subtype counts:\n%s", df_counts)
        logger.info("EC subtype fractions:\n%s", df_frac)
        logger.info("EC subtype fractions (per sample):\n%s", df_frac_sample)

    logger.info("Merge EC zonation subtypes back into vascular set...")
    merged = sub.copy()
    zon = ec.obs["ec_zonation"].astype("object").reindex(merged.obs_names)
    ct = merged.obs[cell_type_col].astype("object")
    merged.obs["cell_type_vascular"] = ct.mask(zon.notna(), "ECs")
    merged.obs["cell_type_incl_zonation"] = zon.fillna(ct.replace("ECs", "otherECs"))
    merged.obs = merged.obs.drop(columns="ec_zonation")

    for col in ["cell_type_vascular", "cell_type_incl_zonation"]:
        merged.obs[col] = pd.Categorical(
            merged.obs[col], categories=list(_constants.cell_type_colors.keys())
        ).remove_unused_categories()
        merged.uns[f"{col}_colors"] = [
            _constants.cell_type_colors[ct] for ct in merged.obs[col].cat.categories
        ]

    plot_embedding(
        merged,
        "X_umap_harmony_10_20",
        [cell_type_col, "cell_type_vascular", "cell_type_incl_zonation"] + zonation[2:],
        plot_path / "vascular_umap.png",
        ncols=6,
    )

    logger.info(f"Exporting spatial plots to: {plot_path / 'spatial'}...")
    if "sample" in merged.obs.columns:
        for s in merged.obs["sample"].unique():
            with plt.rc_context({"figure.figsize": (10, 10)}):
                sc.pl.embedding(
                    merged[merged.obs["sample"] == s],
                    basis="spatial_microns",
                    color="cell_type_incl_zonation",
                    size=9,
                    show=False,
                )
                plt.savefig(plot_path / "spatial" / f"{s}.png", dpi=500)
                plt.close()

    logger.info("Saving to %s", out)
    if "fov" in merged.obs.columns:
        merged.obs["fov"] = merged.obs["fov"].astype(str)

    merged.write(out, compression="gzip")
    logger.info("Done.")


if __name__ == "__main__":
    main()
