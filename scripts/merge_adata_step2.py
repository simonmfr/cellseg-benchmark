import argparse
import json
import logging
import warnings
from os.path import exists
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

from spatialdata import read_zarr

from cellseg_benchmark.adata_utils import (
    dimensionality_reduction,
    dimensionality_reduction_quick,
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

# Load merged (unintegrated) adata
logger.info("Loading adata...")
merged_adata_path = save_path / "adatas" / "adata_merged_temp.h5ad.gz"
if not merged_adata_path.exists():
    logger.error(f"Merged adata not found at: {merged_adata_path}")
    raise FileNotFoundError(f"Merged adata not found at: {merged_adata_path}")
adata = sc.read_h5ad(merged_adata_path)


adata.obsm["spatial"] = adata.obsm.get("spatial_microns", adata.obsm["spatial"])
adata = filter_spatial_outlier_cells(
    adata,
    data_dir=str(base_path),
    sample_paths_file=sample_paths_file,
    save_path=save_path / "plots",
    logger=logger,
)

# Special case: rerun filter_low_quality_cells with lower threshold for vpt_3D
if "vpt_3D" in args.seg_method:
    logger.info("Segmentation method contains 'vpt_3D': applying low-quality cell filtering with min_counts=10 due to smaller cell sizes.")
    adata = filter_low_quality_cells(
        adata,
        save_path=save_path / "plots",
        min_counts=10,
        logger=logger,
    )
else:
    adata = filter_low_quality_cells(
        adata,
        save_path=save_path / "plots",
        logger=logger,
    )
    
adata = filter_genes(adata, save_path=save_path / "plots", logger=logger)

adata = normalize_counts(
    adata, save_path=save_path / "plots", seg_method=args.seg_method, logger=logger
)

# Subset to max_cells
max_cells = 750_000
if adata.n_obs > max_cells:
    logger.info(f"Stratified subsetting from {adata.n_obs:,} to {max_cells:,} cells.")
    rng = np.random.default_rng(seed=42)
    counts = adata.obs["sample"].value_counts()
    frac = max_cells / adata.n_obs
    keep_idx = np.concatenate([
        rng.choice(adata.obs_names[adata.obs["sample"] == s], 
                   size=int(np.floor(c * frac)), replace=False)
        for s, c in counts.items()
    ])
    adata = adata[keep_idx].copy()
    
adata = dimensionality_reduction_quick(adata, save_path=save_path / "plots", logger=logger)
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

logger.info("Deleting temp file...")
final_adata_path = output_path / "adata_integrated.h5ad.gz"
if final_adata_path.exists():
    if temp_adata_path.exists():
        logger.info(f"Integration succeeded. Deleting temp file: {temp_adata_path}")
        temp_adata_path.unlink()
    else:
        logger.warning(f"Temp file not found: {temp_adata_path}")
else:
    logger.error(f"Integration output missing at: {final_adata_path}. Skipping deletion.")

logger.info("Done.")
