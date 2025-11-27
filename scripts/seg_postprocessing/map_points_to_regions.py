import argparse
import logging
from multiprocessing import cpu_count
from os import listdir
from pathlib import Path
from typing import Tuple

from geopandas import read_parquet
from joblib import Parallel, delayed
from pandas import concat
from scanpy import read_h5ad
from tqdm import tqdm

from cellseg_benchmark.adata_utils import plot_spatial_multiplot
from cellseg_benchmark.spatial_mapping import map_points_to_regions_from_anndata

# ---------------------------------------------------------------------
# Worker function for one segmentation method
# ---------------------------------------------------------------------
def process_method(method: str, cohort: str, anatom_annot: dict, data_path: Path) -> None:
    """Process a single segmentation method: read, map, save CSV & plot."""
    method_dir = data_path / "analysis" / cohort / method
    adata_points = read_h5ad(method_dir / "adatas" / "adata_integrated.h5ad.gz")

    # Map points to regions
    results = map_points_to_regions_from_anndata(
        adata_points,
        regions_by_slide=anatom_annot,
        coord_key="spatial_microns",
        slide_key="sample",
        include_boundary=True,
        dissolve=False,
        index_kind="name",
        return_df=True,
    )

    # Combine results into a single DataFrame
    df_all = concat([v["df"] for v in results.values()], ignore_index=True)
    df_all = df_all.set_index("obs_id").reindex(adata_points.obs_names)

    # Add to AnnData
    adata_points.obs["region_mapped"] = df_all["label"].values
    adata_points.obs["region_poly_index"] = df_all["poly_index"].values

    # Save CSV
    df_all.to_csv(method_dir / "spatial_registration.csv")

    # Produce plot
    plot_spatial_multiplot(
        adata_points,
        "region_mapped",
        save_path=method_dir / "plots",
        save_name="spatial_registration.png",
    )

def _process_method_wrapper(method: str, cohort: str, anatom_annot: dict, data_path: Path) -> Tuple[str, str]:
    """Small wrapper so we see failures per method instead of crashing everything."""
    try:
        process_method(method, cohort, anatom_annot, data_path)
        return method, "ok"
    except Exception as e:
        return method, f"error: {e}"

# ---------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------
logger = logging.getLogger("annotation")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
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
args = parser.parse_args()

data_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
if args.seg_methods is not None:
    seg_methods = args.seg_methods
else:
    seg_methods = listdir(data_path / "analysis" / args.cohort)

logger.info("read in reference.")
gdf = read_parquet(
    data_path / "misc" / "brain_regions" / f"{args.cohort}_brain_regions.parquet"
)

anatom_annot = {}
for (sample, label), sub in gdf.groupby(["sample", "label"]):
    anatom_annot.setdefault(sample, {})[label] = list(sub.geometry)

if args.n_jobs == -1:
    n_jobs = cpu_count()
else:
    n_jobs = args.n_jobs

logger.info(
    f"Starting parallel processing of {len(seg_methods)} methods with n_jobs={n_jobs}"
)
results = Parallel(n_jobs=n_jobs)(
    delayed(_process_method_wrapper)(m, args.cohort, anatom_annot, data_path) for m in tqdm(seg_methods, desc="spatial registration")
)
for m, status in results:
    if status != "ok":
        logger.warning(f"method {m}: {status}")
logger.info("Done")
