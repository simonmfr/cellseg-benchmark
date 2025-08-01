import argparse
import json
import logging
import warnings
from pathlib import Path

from spatialdata import read_zarr

from cellseg_benchmark.adata_utils import (
    dimensionality_reduction,
    filter_genes,
    filter_low_quality_cells,
    filter_spatial_outlier_cells,
    integration_harmony,
    merge_adatas,
    normalize_counts,
)

warnings.filterwarnings("ignore")

# Logger setup
logger = logging.getLogger("integrate_adatas")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

# CLI args
parser = argparse.ArgumentParser(
    description="Integrate adatas from a selected segmentation method."
)
parser.add_argument("cohort", help="Cohort name, e.g., 'foxf2'")
parser.add_argument(
    "seg_method", help="Segmentation method, e.g., 'Cellpose_1_nuclei_model'"
)
parser.add_argument(
    "--age",
    default=False,
    action="store_true",
    help="Age information is available. Otherwise, assume age: 6m.",
)
parser.add_argument(
    "--genotype",
    default=False,
    action="store_true",
    help="Genotype information is available. Otherwise, assume WT.",
)
args = parser.parse_args()

# Paths
base_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
samples_path = base_path / "samples"
save_path = base_path / "analysis" / args.cohort / args.seg_method
save_path.mkdir(parents=True, exist_ok=True)

# Load sample paths
with open(base_path / "sample_paths.json") as f:
    sample_paths_file = json.load(f)

excluded_samples = {  # TODO: move to yaml
    "foxf2": ["foxf2_s2_r0", "foxf2_s3_r0", "foxf2_s3_r1"],
    "aging": [],
}.get(args.cohort, [])

# Load sdata
logger.info("Loading data...")
sdata_list = []
available_names = set()

for sample_dir in samples_path.glob(f"{args.cohort}*"):  # restriction to cohort folders
    if sample_dir.name in excluded_samples:
        continue
    sdata = read_zarr(sample_dir / "sdata_z3.zarr", selection=("tables",))
    current_names = ["_".join(k.split("_")[1:]) for k in sdata.tables.keys()]
    available_names.update(current_names)
    sdata_list.append((sample_dir.name, sdata))

# Merge and process
adata = merge_adatas(
    sdata_list,
    seg_method=args.seg_method,
    sample_paths_file=sample_paths_file,
    logger=logger,
    plot_qc_stats=True,
    age=args.age,
    genotype=args.genotype,
    save_path=save_path / "plots",
)
del sdata_list

adata.obsm["spatial"] = adata.obsm.get("spatial_microns", adata.obsm["spatial"])
adata = filter_spatial_outlier_cells(
    adata,
    data_dir=str(base_path),
    sample_paths_file=sample_paths_file,
    save_path=save_path / "plots",
    logger=logger,
)
adata = filter_low_quality_cells(adata, save_path=save_path / "plots", logger=logger)
adata = filter_genes(adata, save_path=save_path / "plots", logger=logger)
adata = normalize_counts(
    adata, save_path=save_path / "plots", seg_method=args.seg_method, logger=logger
)
adata = dimensionality_reduction(adata, save_path=save_path / "plots", logger=logger)
adata = integration_harmony(
    adata, batch_key="sample", save_path=save_path / "plots", logger=logger
)

# Save result
logger.info("Saving integrated object...")
if "fov" not in adata.obs.columns:
    adata.obs["fov"] = ""
adata.obs["fov"] = adata.obs["fov"].astype(str)
output_path = save_path / "adatas"
output_path.mkdir(parents=True, exist_ok=True)
adata.write(output_path / "adata_integrated.h5ad.gz", compression="gzip")
logger.info("Done.")
