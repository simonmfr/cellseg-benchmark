import argparse
import logging
import os
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple

import scanpy as sc
import yaml
from anndata import AnnData
from spatialdata import read_zarr

from cellseg_benchmark.adata_utils import (
    filter_genes,
    filter_low_quality_cells,
    filter_spatial_outlier_cells,
    integration_harmony,
    merge_adatas,
    normalize_counts,
    pca_umap_single,
)


def _load_one(
    sample_dir: Path, seg_method: str, logger: logging.Logger
) -> Tuple[str, AnnData | None]:
    """Load AnnData from one master sdata."""
    sdata = read_zarr(sample_dir / "sdata_z3.zarr", selection=("tables",))
    if f"adata_{seg_method}" not in sdata.tables.keys():
        if logger:
            logger.warning(f"Skipping {seg_method}. No such key: {seg_method}")
        return sample_dir.name, None
    return sample_dir.name, sdata[f"adata_{seg_method}"]


warnings.filterwarnings("ignore", ".*The table is annotating*", UserWarning)
sc.settings.n_jobs = -1

logger = logging.getLogger("integrate_adatas")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

parser = argparse.ArgumentParser(
    description="Integrate adatas from a selected segmentation method."
)
parser.add_argument("cohort", help="Cohort name, e.g., 'foxf2'")
parser.add_argument(
    "seg_method", help="Segmentation method, e.g., 'Cellpose_1_nuclei_model'"
)
args = parser.parse_args()

base_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
samples_path = base_path / "samples"
save_path = base_path / "analysis" / args.cohort / args.seg_method
save_path.mkdir(parents=True, exist_ok=True)

sample_metadata_file, excluded = (
    yaml.safe_load(open(base_path / "misc" / f))
    for f in ["sample_metadata.yaml", "samples_excluded.yaml"]
)
excluded_samples = set(excluded.get(args.cohort, []))

logger.info("Loading data...")
loads = []
for sample_dir in samples_path.glob(f"{args.cohort}*"):  ###############
    if sample_dir.name in excluded_samples:
        continue
    if not (sample_dir / "sdata_z3.zarr").exists():
        logger.error("master sdata in %s not found.", sample_dir)
        continue
    loads.append(sample_dir)

# drop samples where AnnData could not be loaded (missing or invalid), e.g. aging s8 r1 Cellpose 2 Transcripts
with ProcessPoolExecutor(max_workers=int(os.getenv("SLURM_CPUS_PER_TASK", 1))) as ex:
    futures = {ex.submit(_load_one, p, args.seg_method, logger): p for p in loads}
    adata_list = [f.result() for f in as_completed(futures)]
adata_list = [(x, y) for x, y in adata_list if y is not None]

# Merge and process
adata = merge_adatas(
    adata_list,
    seg_method=args.seg_method,
    logger=logger,
    plot_qc_stats=True,
    save_path=save_path / "plots",
)
del adata_list

# workaround to fix adata.obs formatting ###############
adata.obs["sample"] = adata.obs["sample"].str.replace(
    rf"^{args.cohort}_(\d+)_(\d+)$", rf"{args.cohort}_s\1_r\2", regex=True
)
adata.obs["condition"] = (
    adata.obs["genotype"].astype(str) + "_" + adata.obs["age_months"].astype(str)
)
###############

adata.obsm["spatial"] = adata.obsm.get("spatial_microns", adata.obsm["spatial"])
adata = filter_spatial_outlier_cells(
    adata,
    data_dir=str(base_path),
    sample_metadata_file=sample_metadata_file,
    save_path=save_path / "plots",
    logger=logger,
)

# Special case: rerun filter_low_quality_cells with lower threshold for vpt_3D
if "vpt_3D" in args.seg_method:
    logger.info(
        "Segmentation method contains 'vpt_3D': applying low-quality cell filtering with min_counts=10 due to smaller cell sizes."
    )
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

adata = pca_umap_single(adata, save_path=save_path / "plots", logger=logger)

adata = integration_harmony(
    adata,
    batch_key="slide",
    save_path=save_path / "plots",
    logger=logger,
)

logger.info("Saving integrated object...")
if "fov" not in adata.obs.columns:
    adata.obs["fov"] = ""
adata.obs["fov"] = adata.obs["fov"].astype(str)
output_path = save_path / "adatas"
output_path.mkdir(parents=True, exist_ok=True)
adata.write(output_path / "adata_integrated.h5ad.gz", compression="gzip")

logger.info("Done.")
