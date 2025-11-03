import argparse
import logging
from os import listdir
from pathlib import Path

from geopandas import read_parquet
from pandas import concat
from scanpy import read_h5ad
from tqdm import tqdm

from cellseg_benchmark.spatial_mapping import map_points_to_regions_from_anndata
from cellseg_benchmark.adata_utils import plot_spatial_multiplot

logger = logging.getLogger("annotation")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

parser = argparse.ArgumentParser(
    description="Map points to brain regions."
)
parser.add_argument("cohort", help="cohort name")
parser.add_argument("seg_methods", nargs=argparse.REMAINDER,
                    help="If not set, all segmentation methods within cohort will be mapped. Otherwise provide segmentation method to annotate."
                    )
args = parser.parse_args()

data_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
if args.seg_methods is not None:
    seg_methods = args.seg_methods
else:
    seg_methods = listdir(data_path / "analysis" / args.cohort)

logger.info("read in reference.")
gdf = read_parquet(data_path / "misc" / "brain_regions" / f"{args.cohort}_brain_regions.parquet")
anatom_annot = {}
for (sample, label), sub in gdf.groupby(["sample", "label"]):
    anatom_annot.setdefault(sample, {})[label] = list(sub.geometry)

for method in tqdm(seg_methods):
    logger.info(f"processing {method}")
    adata_points = read_h5ad(data_path / "analysis" / args.cohort / method / "adatas" / "adata_integrated.h5ad.gz")
    results = map_points_to_regions_from_anndata(adata_points,
                                                 regions_by_slide=anatom_annot,
                                                 coord_key="spatial_microns",
                                                 slide_key="sample",
                                                 include_boundary=True,
                                                 dissolve=False,
                                                 index_kind="name",
                                                 return_df=True
                                                 )
    df_all = concat([v['df'] for v in results.values()], ignore_index=True)
    df_all = df_all.set_index("obs_id").reindex(adata_points.obs_names)
    adata_points.obs["region_mapped"] = df_all['label'].values
    adata_points.obs["region_poly_index"] = df_all['poly_index'].values
    df_all.to_csv(data_path / "analysis" / args.cohort / method / "spatial_registration.csv")
    plot_spatial_multiplot(adata_points, "region_mapped",
                           save_path=data_path / "analysis" / args.cohort / method,
                           save_name="spatial_registration.png",
                           sort=True
                           )

