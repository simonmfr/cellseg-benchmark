#!/usr/bin/env python
"""Turn joint BANKSY clusters into anatomical brain-region polygons.

The clustering is joint across all samples of a cohort, so cluster k means the
same thing in every section. The cluster -> region assignment is therefore one
global table of ~20 entries in a YAML config, not one decision per polygon per
sample. Per-sample exceptions stay possible but have to be written down.

Two steps:

1. ``--init`` writes a YAML skeleton plus, per cluster, a multi-sample plot and
   a marker-expression table, i.e. the evidence needed to name the clusters.
2. Without ``--init`` the filled-in YAML is applied and the region parquet
   consumed by ``map_points_to_regions.py`` is written.
"""

import argparse
import logging
import pathlib

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from shapely.geometry import Point
from skimage.measure import label as cc_label

from cellseg_benchmark._constants import brain_regions_colors
from cellseg_benchmark.adata_utils import plot_spatial_multiplot
from cellseg_benchmark.spatial_mapping import (
    _extent_from_geo,
    gridify,
    polygons_per_component_exact,
    relabel_small_holes,
    remove_small_islands,
    to_gdf,
)

BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

# Shown next to each cluster in --init so clusters are named from evidence
# rather than from shape alone. Genes absent from the panel are skipped.
MARKERS = {
    "DG-sg": ["Prox1", "Dock10"],
    "CAsp": ["Fibcd1", "Wfs1", "Neurod6"],
    "CTX": ["Cux2", "Rorb", "Foxp2", "Bcl11b"],
    "STR": ["Ppp1r1b", "Drd1", "Adora2a"],
    "BS": ["Tcf7l2", "Slc17a6"],
    "fiber_tracts": ["Mbp", "Plp1", "Mog"],
    "VS": ["Ttr", "Foxj1"],
    "Meninges": ["Dcn", "Slc47a1"],
}

DEFAULT_CLEANUP = {"min_hole_area_um2": 150000.0, "min_island_area_um2": 60000.0}

logger = logging.getLogger("brain_regions")
logger.setLevel(logging.INFO)
_handler = logging.StreamHandler()
_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(_handler)


def cluster_codes(adata, cluster_key):
    """Factorize the cluster column once, so codes mean the same in every sample.

    Args:
        adata: AnnData holding all samples of the cohort.
        cluster_key: `.obs` column with the BANKSY clusters.

    Returns:
        Mapping from integer code to cluster id. The codes are stored in
        `.obs["_cluster_code"]`.
    """
    codes, uniques = pd.factorize(adata.obs[cluster_key].astype(str), sort=True)
    adata.obs["_cluster_code"] = codes
    return {i: str(u) for i, u in enumerate(uniques)}


def write_skeleton(adata, cluster_key, config_path, plot_dir):
    """Write a YAML skeleton and the per-cluster evidence needed to fill it in.

    Args:
        adata: AnnData holding all samples of the cohort.
        cluster_key: `.obs` column with the BANKSY clusters.
        config_path: Where to write the YAML skeleton.
        plot_dir: Folder for the per-cluster plots and the marker table.
    """
    plot_dir.mkdir(parents=True, exist_ok=True)
    clusters = sorted(adata.obs[cluster_key].astype(str).unique(), key=_sort_key)

    genes = [g for gs in MARKERS.values() for g in gs if g in adata.var_names]
    if genes:
        expr = sc.get.obs_df(
            adata, keys=genes + [cluster_key], layer="volume_log1p_norm"
        )
        table = expr.groupby(cluster_key, observed=True).mean().round(3)
        table.to_csv(plot_dir / "cluster_markers.csv")
        logger.info(
            "Marker means for %d/%d genes -> cluster_markers.csv",
            len(genes),
            sum(map(len, MARKERS.values())),
        )
    else:
        logger.warning("No marker genes of MARKERS are in the panel.")

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
    skeleton = {
        "cluster_key": cluster_key,
        "clusters": {c: None for c in clusters},
        "sample_overrides": {},
        "point_overrides": {},
        "cleanup": {"default": dict(DEFAULT_CLEANUP)},
    }
    with open(config_path, "w") as fh:
        yaml.safe_dump(skeleton, fh, sort_keys=False, default_flow_style=False)
    logger.info(
        "Wrote %s: fill in one region per cluster, then rerun without --init.",
        config_path,
    )


def _sort_key(c):
    """Sort cluster ids numerically when possible.

    Args:
        c: Cluster id.

    Returns:
        A tuple ordering numeric ids before non-numeric ones.
    """
    return (0, int(c)) if str(c).isdigit() else (1, str(c))


def build_regions(adata, cfg, code_to_cluster, plot_dir):
    """Clean the cluster raster per sample and label its connected components.

    Args:
        adata: AnnData holding all samples of the cohort.
        cfg: Parsed YAML config.
        code_to_cluster: Mapping from integer code to cluster id.
        plot_dir: Folder for the per-sample QC figure.

    Returns:
        `{sample: {label: [Polygon, ...]}}`.
    """
    mapping = {str(k): v for k, v in (cfg.get("clusters") or {}).items()}
    unnamed = [c for c, v in mapping.items() if not v]
    if unnamed:
        logger.warning("Clusters without a region, dropped: %s", ", ".join(unnamed))

    cleanup = cfg.get("cleanup") or {}
    default_cleanup = {**DEFAULT_CLEANUP, **(cleanup.get("default") or {})}

    out, grids = {}, {}
    for sample in sorted(adata.obs["sample"].astype(str).unique()):
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

        # cluster -> region label, with per-sample exceptions.
        overrides = {
            str(k): v
            for k, v in ((cfg.get("sample_overrides") or {}).get(sample) or {}).items()
        }
        for reg in regions:
            reg["label"] = overrides.get(reg["value"], mapping.get(reg["value"]))

        # A cluster that is genuinely two anatomies inside one section is split by
        # naming a point in the offending component, which survives re-runs.
        for x, y, label in (cfg.get("point_overrides") or {}).get(sample, []):
            hits = [r for r in regions if r["poly"].covers(Point(float(x), float(y)))]
            if not hits:
                logger.warning(
                    "%s: point override (%s, %s) hits no region.", sample, x, y
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
    """Render every sample's labelled raster as one multi-panel QC figure.

    Args:
        grids: Mapping from sample to `(clean_grid, geo, regions)`.
        plot_dir: Folder to write the figure into.
        n_cols: Panels per row.
    """
    labels = sorted(
        {r["label"] for _, _, regs in grids.values() for r in regs if r["label"]}
    )
    colors = [_color(lab, i) for i, lab in enumerate(labels)]
    lut = {lab: i for i, lab in enumerate(labels)}

    n_rows = int(np.ceil(len(grids) / n_cols))
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
    for ax, (sample, (clean, geo, regions)) in zip(np.ravel(axs), grids.items()):
        # Colour per connected component, not per cluster code: a point override
        # relabels one component only.
        label_of = {(r["code"], r["comp_id"]): r["label"] for r in regions}
        img = np.full(clean.shape, np.nan)
        for code in np.unique(clean[clean >= 0]):
            cc = cc_label(clean == code, connectivity=1)
            for cid in range(1, cc.max() + 1):
                label = label_of.get((int(code), cid - 1))
                if label:
                    img[cc == cid] = lut[label]
        ax.imshow(
            img,
            origin="upper",
            extent=_extent_from_geo(clean, geo),
            cmap=ListedColormap(colors),
            vmin=-0.5,
            vmax=len(labels) - 0.5,
            interpolation="nearest",
        )
        ax.set_title(sample)
        ax.set_aspect("equal")
        ax.axis("off")
    for ax in np.ravel(axs)[len(grids) :]:
        ax.axis("off")

    fig.legend(
        handles=[Patch(facecolor=c, label=lab) for lab, c in zip(labels, colors)],
        loc="center right",
    )
    plot_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_dir / "brain_regions.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def _color(label, i):
    """Look up a region colour, falling back to tab20.

    Args:
        label: Region label.
        i: Position of the label, used for the fallback colour.

    Returns:
        A matplotlib colour.
    """
    return brain_regions_colors.get(label, plt.get_cmap("tab20")(i % 20))


def main():
    """Parse arguments and either write the YAML skeleton or the region parquet."""
    parser = argparse.ArgumentParser(description=__doc__)
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
    plot_dir = pathlib.Path(args.adata_path).parent / "plots" / "brain_regions"

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

    code_to_cluster = cluster_codes(adata, cluster_key)
    regions_by_slide = build_regions(adata, cfg, code_to_cluster, plot_dir)

    out = pathlib.Path(
        args.out
        or BASE_PATH / "misc" / "brain_regions" / f"{args.cohort}_brain_regions.parquet"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    to_gdf(regions_by_slide).to_parquet(out)
    logger.info("Wrote %s", out)


if __name__ == "__main__":
    main()
