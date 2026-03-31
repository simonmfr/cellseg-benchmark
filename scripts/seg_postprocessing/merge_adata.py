import argparse
import logging
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from typing import Tuple

import pandas as pd
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
yaml_samples = [
    name
    for name, meta in sample_metadata_file.items()
    if meta.get("cohort") == args.cohort and name not in excluded_samples
]
if (
    args.cohort == "htra1"
):  # add 6/18m WT samples from aging cohort as additional controls
    yaml_samples += [
        "aging_s1_r1",
        "aging_s5_r1",
        "aging_s6_r0",
        "aging_s7_r1",
        "aging_s7_r2",
        "aging_s8_r2",
        "aging_s11_r0",
    ]

logger.info("Loading data...")
loads = []
for name in yaml_samples:
    p = samples_path / name
    if not (p / "sdata_z3.zarr").exists():
        logger.error("master sdata in %s not found.", p)
        continue
    loads.append(p)

max_workers = int(os.getenv("SLURM_CPUS_PER_TASK", 1))
loader = partial(_load_one, seg_method=args.seg_method, logger=None)

with ProcessPoolExecutor(max_workers=max_workers) as ex:
    results = list(ex.map(loader, loads))

# keep YAML order, drop Nones
adata_list = [(name, adata) for name, adata in results if adata is not None]

# temp fix for aging_s11_r0
for i, (n, ad) in enumerate(adata_list):
    if n == "aging_s11_r0":
        # if "aging" in n:
        m = sample_metadata_file[n]
        for k, v in {
            **m,
            "sample": n,
            "condition": f"{m['genotype']}_{m['age_months']}",
        }.items():
            ad.obs[k] = pd.Categorical([str(v)] * len(ad))
        adata_list[i] = (n, ad)

# Merge and process
adata = merge_adatas(
    adata_list,
    seg_method=args.seg_method,
    logger=logger,
    plot_qc_stats=True,
    save_path=save_path / "plots",
)
del adata_list

adata.obsm["spatial"] = adata.obsm.get("spatial_microns", adata.obsm["spatial"])
if args.seg_method == "Cellpose_1_Merlin":  # workaround, as explorer is in pixel units
    adata.obsm["spatial"] = adata.obsm.get("spatial_pixel", adata.obsm["spatial"])

adata = filter_spatial_outlier_cells(
    adata,
    data_dir=str(base_path),
    sample_metadata_file=sample_metadata_file,
    save_path=save_path / "plots",
    logger=logger,
)

if args.seg_method == "Cellpose_1_Merlin":  # workaround, as explorer is in pixel units
    adata.obsm["spatial"] = adata.obsm.get("spatial_microns", adata.obsm["spatial"])

if "vpt_3D" in args.seg_method:  # min_counts=10 due to smaller cell sizes
    min_counts = 10
elif args.cohort == "SynergyLung":  # more lenient for initial analysis
    min_counts = 15
else:
    min_counts = None  # default = 25

adata = filter_low_quality_cells(
    adata,
    save_path=save_path / "plots",
    **({"min_counts": min_counts} if min_counts is not None else {}),
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
