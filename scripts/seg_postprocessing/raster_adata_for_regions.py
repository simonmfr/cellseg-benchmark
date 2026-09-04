#!/usr/bin/env python
"""Merge the raster segmentation of one cohort into a reference for region mapping.

Unlike merge_adata.py this keeps every bin lying on tissue. Cell-level QC drops
low-count bins, which in a 25 um raster means white matter, ventricles and
meninges, and the resulting holes would be inherited by every segmentation
method mapped against these regions. Only anatomically shifted bins are removed.
"""

import argparse
import logging
import pathlib
import warnings

import spatialdata as sd
import yaml

from cellseg_benchmark import _constants, adata_utils

warnings.filterwarnings("ignore", ".*The table is annotating*", UserWarning)
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s]: %(message)s"
)
logger = logging.getLogger("raster_regions")

parser = argparse.ArgumentParser(description="Raster adata for brain-region mapping.")
parser.add_argument("cohort", help="Cohort name, e.g. 'aging'.")
parser.add_argument("--seg_method", default="Negative_Control_Rastered_25")
args = parser.parse_args()

base_path = pathlib.Path(_constants.BASE_PATH)
save_path = base_path / "analysis" / args.cohort / args.seg_method
(save_path / "plots").mkdir(parents=True, exist_ok=True)

sample_metadata_file, excluded = (
    yaml.safe_load(open(base_path / "misc" / f))
    for f in ["sample_metadata.yaml", "samples_excluded.yaml"]
)
excluded_samples = set(excluded.get(args.cohort, []))
samples = [
    name
    for name, meta in sample_metadata_file.items()
    if meta.get("cohort") == args.cohort and name not in excluded_samples
]
if args.cohort == "htra1":
    samples += _constants.htra1_aging_controls

adata_list = []
for name in samples:
    path = base_path / "samples" / name / "sdata_z3.zarr"
    if not path.exists():
        logger.error("%s: no master sdata", name)
        continue
    sdata = sd.read_zarr(path, selection=("tables",))
    key = f"adata_{args.seg_method}"
    if key in sdata.tables:
        adata_list.append((name, sdata[key]))
    else:
        logger.warning("%s: %s missing", name, key)

logger.info("Merging %d samples...", len(adata_list))
adata = adata_utils.merge_adatas(
    adata_list, seg_method=args.seg_method, logger=logger, plot_qc_stats=False
)
del adata_list
adata.obsm["spatial"] = adata.obsm["spatial_microns"]

adata = adata_utils.filter_spatial_outlier_cells(
    adata,
    data_dir=str(base_path),
    sample_metadata_file=sample_metadata_file,
    save_path=save_path / "plots",
    logger=logger,
)
adata = adata_utils.filter_genes(adata, save_path=save_path / "plots", logger=logger)
adata = adata_utils.normalize_counts(
    adata,
    save_path=save_path / "plots",
    seg_method=args.seg_method,
    trim_outliers=False,
    logger=logger,
)

out = save_path / "adatas" / "adata_regions.h5ad.gz"
out.parent.mkdir(parents=True, exist_ok=True)
adata.write(out, compression="gzip")
logger.info(
    "Wrote %s: %d bins, %d genes, %d samples",
    out,
    adata.n_obs,
    adata.n_vars,
    adata.obs["sample"].nunique(),
)
