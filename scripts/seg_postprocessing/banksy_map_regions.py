#!/usr/bin/env python
import argparse
import logging
import os
import pathlib
import sys
import typing

import geopandas as gpd
import pandas as pd
import scanpy as sc
import tqdm
from joblib import Parallel, delayed

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))

from cellseg_benchmark import BASE_PATH
from cellseg_benchmark._constants import brain_regions_broad, brain_regions_colors
from cellseg_benchmark.adata_utils import plot_spatial_multiplot
from cellseg_benchmark.spatial_mapping import map_points_to_regions_from_anndata


def process_method(
    method: str, cohort: str, anatom_annot: dict, data_path: pathlib.Path, title_keys
) -> None:
    """Process a single segmentation method: read, map, save CSV & plot."""
    method_dir = data_path / "analysis" / cohort / method
    adata_points = sc.read_h5ad(method_dir / "adatas" / "adata_integrated.h5ad.gz")

    results = map_points_to_regions_from_anndata(
        adata_points,
        regions_by_slide=anatom_annot,
        coord_key="spatial_microns",
        slide_key="sample",
        include_boundary=True,
        index_kind="name",
        return_df=True,
    )

    df_all = pd.concat([v["df"] for v in results.values()], ignore_index=True)
    df_all = df_all.set_index("obs_id").reindex(adata_points.obs_names)
    df_all = df_all.rename(columns={"label": "brain_region"})
    df_all["brain_region_broad"] = df_all["brain_region"].map(
        lambda lab: brain_regions_broad.get(lab, lab)
    )

    adata_points.obs["brain_region"] = df_all["brain_region"].values
    adata_points.obs["brain_region_broad"] = df_all["brain_region_broad"].values
    adata_points.obs["brain_region_poly_index"] = df_all["poly_index"].values

    df_all.to_csv(method_dir / "brain_regions.csv")

    plot_spatial_multiplot(
        adata_points,
        "brain_region",
        save_path=method_dir / "plots",
        save_name="brain_regions.png",
        palette=brain_regions_colors,
        title_keys=title_keys,
    )


def _process_method_wrapper(
    method: str, cohort: str, anatom_annot: dict, data_path: pathlib.Path, title_keys
) -> typing.Tuple[str, str]:
    """Small wrapper so we see failures per method instead of crashing everything."""
    try:
        process_method(method, cohort, anatom_annot, data_path, title_keys)
        return method, "ok"
    except Exception as e:
        return method, f"error: {e}"


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s]: %(message)s"
)
logger = logging.getLogger("annotation")

parser = argparse.ArgumentParser(description="Map points to brain regions.")
parser.add_argument("cohort", help="cohort name")
parser.add_argument(
    "--seg_methods",
    nargs=argparse.REMAINDER,
    help=(
        "If not set, all segmentation methods within cohort will be mapped. "
        "Otherwise provide segmentation methods to annotate."
    ),
)
parser.add_argument(
    "--n_jobs",
    type=int,
    default=-1,
    help="Number of parallel jobs for segmentation methods (-1 = all cores).",
)
parser.add_argument(
    "--title_keys",
    nargs="+",
    default=["sample"],
    help="obs column(s) for panel titles, e.g. --title_keys sample age.",
)
args = parser.parse_args()

data_path = pathlib.Path(BASE_PATH)
if args.seg_methods is not None:
    seg_methods = args.seg_methods
else:
    seg_methods = os.listdir(data_path / "analysis" / args.cohort)

logger.info("read in reference.")
gdf = gpd.read_parquet(
    data_path / "misc" / "brain_regions" / f"{args.cohort}_brain_regions.parquet"
)

anatom_annot = {}
for (sample, label), sub in gdf.groupby(["sample", "label"]):
    anatom_annot.setdefault(sample, {})[label] = list(sub.geometry)

logger.info(
    f"Starting parallel processing of {len(seg_methods)} methods with "
    f"n_jobs={args.n_jobs}"
)
results = Parallel(n_jobs=args.n_jobs)(
    delayed(_process_method_wrapper)(
        m, args.cohort, anatom_annot, data_path, args.title_keys
    )
    for m in tqdm.tqdm(seg_methods, desc="brain regions")
)
for m, status in results:
    if status != "ok":
        logger.warning(f"method {m}: {status}")
logger.info("Done")
