#!/usr/bin/env python
"""BANKSY joint clustering -> anatomical brain-region polygons -> per-cell labels.

Full cohort pipeline, in order:
    1. banksy_prep_raster.py {cohort}                 -> reference adata (rastered 25um bins)
    2. banksy_clustering.py {cohort} <adata> <dir>    -> joint BANKSY clusters
       (steps 1+2 together: sbatch_utils/banksy_script_creation.py {cohort})
    3. this script --init                             -> YAML skeleton + per-cluster evidence
    4. fill in configs/brain_regions/{cohort}.yaml, rerun this script without --init
    5. banksy_map_regions.py {cohort}                 -> per-cell labels, one method's
       adatas/adata_integrated.h5ad.gz at a time

Clusters are shared across the whole cohort, so one YAML table names every
cluster once. Per-sample exceptions go in sample_overrides/point_overrides
instead.

Writes, per run: {plot_dir}/brain_regions.png (final render), components.png
+ components.csv (QC: per-component id/coords for point_overrides), and the
region parquet with both a fine `label` and coarse `label_broad` column
(cellseg_benchmark._constants.brain_regions_broad).

    # 1. YAML skeleton + per-cluster plots and marker table to name them from
    banksy_annotate_regions.py aging adata_regions.h5ad.gz --cluster-key banksy_coarse_k50_res0.4 --init
    # 2. apply the filled-in YAML: cleans holes/islands, writes the region parquet
    banksy_annotate_regions.py aging adata_regions.h5ad.gz
"""

import argparse
import logging
import pathlib
import re
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from shapely.geometry import Point
from skimage.measure import label as cc_label
from skimage.segmentation import find_boundaries

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))

from cellseg_benchmark import BASE_PATH
from cellseg_benchmark._constants import brain_regions_broad, brain_regions_colors
from cellseg_benchmark._markers import brain_region_markers
from cellseg_benchmark.adata_utils import plot_spatial_multiplot
from cellseg_benchmark.spatial_mapping import (
    _extent_from_geo,
    gridify,
    polygons_per_component_exact,
    relabel_small_holes,
    remove_small_islands,
    to_gdf,
)

BASE_PATH = pathlib.Path(BASE_PATH)
DEFAULT_CLEANUP = {"min_hole_area_um2": 300000.0, "min_island_area_um2": 30000.0}

logger = logging.getLogger("brain_regions")
logger.setLevel(logging.INFO)
_handler = logging.StreamHandler()
_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(_handler)


def write_skeleton(adata, cluster_key, config_path, plot_dir):
    """Write a YAML skeleton plus the per-cluster evidence needed to fill it in."""
    plot_dir.mkdir(parents=True, exist_ok=True)
    clusters = sorted(adata.obs[cluster_key].astype(str).unique(), key=int)

    genes = [
        g for gs in brain_region_markers.values() for g in gs if g in adata.var_names
    ]
    if genes:
        expr = sc.get.obs_df(
            adata, keys=genes + [cluster_key], layer="volume_log1p_norm"
        )
        expr.groupby(cluster_key, observed=True).mean().round(3).to_csv(
            plot_dir / "cluster_markers.csv"
        )
        logger.info("%d marker genes in the panel -> cluster_markers.csv", len(genes))
    else:
        logger.warning("None of the marker genes are in the panel.")

    for c in clusters:
        adata.obs["_one"] = np.where(
            adata.obs[cluster_key].astype(str) == c, c, "other"
        )
        plot_spatial_multiplot(
            adata,
            "_one",
            save_path=str(plot_dir),
            save_name=f"cluster_{c}.png",
            palette={c: "#d62728", "other": "#e0e0e0"},
            title=f"{cluster_key} = {c}",
        )
    del adata.obs["_one"]

    config_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_path, "w") as fh:
        yaml.safe_dump(
            {
                "cluster_key": cluster_key,
                "clusters": {c: None for c in clusters},
                "sample_overrides": {},
                "point_overrides": {},
                "cleanup": {"default": dict(DEFAULT_CLEANUP)},
            },
            fh,
            sort_keys=False,
        )
    logger.info("Wrote %s: name every cluster, then rerun without --init.", config_path)


def build_regions(adata, cfg, code_to_cluster, plot_dir):
    """Clean the cluster raster per sample and label its connected components."""
    mapping = {str(k): v for k, v in (cfg.get("clusters") or {}).items()}
    unnamed = [c for c, v in mapping.items() if not v]
    if unnamed:
        logger.warning("Clusters without a region, dropped: %s", ", ".join(unnamed))

    cleanup = cfg.get("cleanup") or {}
    default_cleanup = {**DEFAULT_CLEANUP, **(cleanup.get("default") or {})}

    out, grids = {}, {}
    samples = sorted(
        adata.obs["sample"].astype(str).unique(),
        key=lambda s: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)],
    )
    for sample in samples:
        sub = adata[adata.obs["sample"].astype(str) == sample]
        grid, geo = gridify(sub, "_cluster_code")
        _, _, dx, dy = geo
        thresholds = {**default_cleanup, **(cleanup.get(sample) or {})}

        clean = relabel_small_holes(
            grid, min_hole_area_um2=thresholds["min_hole_area_um2"], dx=dx, dy=dy
        )
        if thresholds.get("min_island_area_um2"):
            clean = remove_small_islands(
                clean,
                min_island_area_um2=thresholds["min_island_area_um2"],
                dx=dx,
                dy=dy,
            )

        regions = polygons_per_component_exact(clean, geo, value_map=code_to_cluster)

        overrides = {
            str(k): v
            for k, v in ((cfg.get("sample_overrides") or {}).get(sample) or {}).items()
        }
        for reg in regions:
            reg["label"] = overrides.get(reg["value"], mapping.get(reg["value"]))

        # A cluster that is two anatomies inside one section is split by naming a
        # point in the offending component; unlike an index this survives reruns.
        for x, y, label in (cfg.get("point_overrides") or {}).get(sample, []):
            hits = [r for r in regions if r["poly"].covers(Point(float(x), float(y)))]
            if not hits:
                logger.warning(
                    "%s: point override (%s, %s) hits nothing.", sample, x, y
                )
            for r in hits:
                r["label"] = label

        out[sample] = {}
        for reg in regions:
            if reg["label"]:
                out[sample].setdefault(reg["label"], []).append(reg["poly"])
        grids[sample] = (clean, geo, regions)
        logger.info(
            "%s: %d components -> %d regions", sample, len(regions), len(out[sample])
        )

    _plot_region_grids(grids, plot_dir)
    return out


def _plot_region_grids(grids, plot_dir, n_cols=3):
    """Render every sample's labelled raster as one multi-panel QC figure."""
    present = {r["label"] for _, _, regs in grids.values() for r in regs if r["label"]}
    labels = [lab for lab in brain_regions_colors if lab in present]
    labels += sorted(present.difference(labels))
    cmap = plt.get_cmap("tab20")
    colors = [
        brain_regions_colors.get(lab, cmap(i % 20)) for i, lab in enumerate(labels)
    ]
    lut = {lab: i for i, lab in enumerate(labels)}

    n_rows = int(np.ceil(len(grids) / n_cols))
    fig_final, axs_final = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
    fig_qc, axs_qc = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
    coord_rows = []
    for ax_final, ax_qc, (sample, (clean, geo, regions)) in zip(
        np.ravel(axs_final), np.ravel(axs_qc), grids.items()
    ):
        label_of = {(r["code"], r["comp_id"]): r["label"] for r in regions}
        img = np.full(clean.shape, np.nan)
        comp_ids = np.full(clean.shape, -1, dtype=int)
        next_id = 0
        for code in np.unique(clean[clean >= 0]):
            cc = cc_label(clean == code, connectivity=1)
            for cid in range(1, cc.max() + 1):
                label = label_of.get((int(code), cid - 1))
                if label:
                    img[cc == cid] = lut[label]
                    comp_ids[cc == cid] = next_id
                    next_id += 1

        extent = _extent_from_geo(clean, geo)
        for ax in (ax_final, ax_qc):
            ax.imshow(
                img,
                origin="upper",
                extent=extent,
                cmap=ListedColormap(colors),
                vmin=-0.5,
                vmax=len(labels) - 0.5,
                interpolation="nearest",
            )
            ax.set_title(sample)
            ax.set_aspect("equal")
            ax.axis("off")

        boundary = find_boundaries(comp_ids, mode="outer")
        overlay = np.zeros((*clean.shape, 4))
        overlay[boundary] = (0, 0, 0, 0.5)
        ax_qc.imshow(overlay, origin="upper", extent=extent)
        for i, reg in enumerate(regions):
            c = reg["poly"].representative_point()
            if reg["label"]:
                ax_qc.annotate(
                    str(i),
                    (c.x, c.y),
                    ha="center",
                    va="center",
                    fontsize=6,
                    bbox=dict(
                        boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7
                    ),
                )
            coord_rows.append(
                {
                    "sample": sample,
                    "id": i,
                    "x": round(c.x),
                    "y": round(c.y),
                    "label": reg["label"],
                }
            )
    for ax in [*np.ravel(axs_final)[len(grids) :], *np.ravel(axs_qc)[len(grids) :]]:
        ax.axis("off")

    legend_handles = [Patch(facecolor=c, label=lab) for lab, c in zip(labels, colors)]
    fig_final.legend(handles=legend_handles, loc="center right", title="Brain region")
    fig_qc.legend(handles=legend_handles, loc="center right", title="Brain region")

    plot_dir.mkdir(parents=True, exist_ok=True)
    fig_qc.savefig(plot_dir / "components.png", dpi=150, bbox_inches="tight")
    plt.close(fig_qc)
    pd.DataFrame(coord_rows).to_csv(plot_dir / "components.csv", index=False)
    fig_final.savefig(plot_dir / "brain_regions.png", dpi=150, bbox_inches="tight")
    plt.close(fig_final)


def main():
    """Write the YAML skeleton, or apply it and write the region parquet."""
    parser = argparse.ArgumentParser(description="BANKSY clusters -> brain regions.")
    parser.add_argument("cohort", help="Cohort name, e.g. 'aging'.")
    parser.add_argument("adata_path", help="Adata written by banksy_clustering.py.")
    parser.add_argument(
        "--cluster-key", help="BANKSY cluster column. Required with --init."
    )
    parser.add_argument(
        "--config", help="YAML config. Default: configs/brain_regions/{cohort}.yaml"
    )
    parser.add_argument(
        "--out",
        help="Output parquet. Default: misc/brain_regions/{cohort}_brain_regions.parquet",
    )
    parser.add_argument(
        "--init",
        action="store_true",
        help="Write a YAML skeleton and per-cluster evidence instead of the parquet.",
    )
    args = parser.parse_args()

    repo = pathlib.Path(__file__).resolve().parents[2]
    config = pathlib.Path(
        args.config or repo / "configs" / "brain_regions" / f"{args.cohort}.yaml"
    )
    plot_dir = pathlib.Path(args.adata_path).parent / "plots"

    logger.info("Loading %s", args.adata_path)
    adata = sc.read_h5ad(args.adata_path)

    if args.init:
        if not args.cluster_key:
            parser.error("--init requires --cluster-key")
        write_skeleton(adata, args.cluster_key, config, plot_dir)
        return

    with open(config) as fh:
        cfg = yaml.safe_load(fh)
    cluster_key = args.cluster_key or cfg["cluster_key"]
    logger.info("Using %s from %s", cluster_key, config)

    # Factorize once, so a code means the same cluster in every sample.
    codes, uniques = pd.factorize(adata.obs[cluster_key].astype(str), sort=True)
    adata.obs["_cluster_code"] = codes
    code_to_cluster = dict(enumerate(map(str, uniques)))

    regions_by_slide = build_regions(adata, cfg, code_to_cluster, plot_dir)

    out = pathlib.Path(
        args.out
        or BASE_PATH / "misc" / "brain_regions" / f"{args.cohort}_brain_regions.parquet"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    gdf = to_gdf(regions_by_slide)
    gdf["label_broad"] = gdf["label"].map(lambda lab: brain_regions_broad.get(lab, lab))
    gdf.to_parquet(out)
    logger.info("Wrote %s", out)


if __name__ == "__main__":
    main()
