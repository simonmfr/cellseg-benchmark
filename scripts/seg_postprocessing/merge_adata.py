import argparse
import concurrent.futures
import functools
import logging
import os
import typing
import warnings
from pathlib import Path

import anndata
import scanpy as sc
import spatialdata as sd
import yaml

from cellseg_benchmark import _constants, adata_utils

warnings.filterwarnings("ignore", ".*The table is annotating*", UserWarning)
sc.settings.n_jobs = -1

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s]: %(message)s"
)
logger = logging.getLogger("integrate_adatas")


def _load_one(
    sample_dir: Path, seg_method: str, logger: logging.Logger
) -> typing.Tuple[str, anndata.AnnData | None]:
    """Load AnnData from one master sdata."""
    sdata = sd.read_zarr(sample_dir / "sdata_z3.zarr", selection=("tables",))
    if f"adata_{seg_method}" not in sdata.tables.keys():
        if logger:
            logger.warning(f"Skipping {seg_method}. No such key: {seg_method}")
        return sample_dir.name, None
    return sample_dir.name, sdata[f"adata_{seg_method}"]


def main():
    """Integrate adatas from a selected segmentation method."""
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("cohort", help="Cohort name, e.g., 'foxf2'")
    parser.add_argument(
        "seg_method", help="Segmentation method, e.g., 'Cellpose_1_nuclei_model'"
    )
    args = parser.parse_args()

    base_path = Path(_constants.BASE_PATH)
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
    ):  # add WT samples from aging cohort as additional controls
        yaml_samples += _constants.htra1_aging_controls

    logger.info("Loading data...")
    loads = []
    for name in yaml_samples:
        p = samples_path / name
        if not (p / "sdata_z3.zarr").exists():
            logger.error("master sdata in %s not found.", p)
            continue
        loads.append(p)

    max_workers = int(os.getenv("SLURM_CPUS_PER_TASK", 1))
    loader = functools.partial(_load_one, seg_method=args.seg_method, logger=None)

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
        results = list(ex.map(loader, loads))

    # keep YAML order, drop Nones
    adata_list = [(name, ad) for name, ad in results if ad is not None]

    # Merge and process
    adata = adata_utils.merge_adatas(
        adata_list,
        seg_method=args.seg_method,
        logger=logger,
        plot_qc_stats=True,
        save_path=save_path / "plots",
    )
    del adata_list

    adata.obsm["spatial"] = adata.obsm.get("spatial_microns", adata.obsm["spatial"])
    if (
        args.seg_method == "Cellpose_1_Merlin"
    ):  # workaround, as explorer is in pixel units
        adata.obsm["spatial"] = adata.obsm.get("spatial_pixel", adata.obsm["spatial"])

    adata = adata_utils.filter_spatial_outlier_cells(
        adata,
        data_dir=str(base_path),
        sample_metadata_file=sample_metadata_file,
        save_path=save_path / "plots",
        logger=logger,
    )

    if (
        args.seg_method == "Cellpose_1_Merlin"
    ):  # workaround, as explorer is in pixel units
        adata.obsm["spatial"] = adata.obsm.get("spatial_microns", adata.obsm["spatial"])

    if "vpt_3D" in args.seg_method:
        min_counts = 10  # min_counts = 10 due to smaller cell sizes
    else:
        min_counts = None  # default = 25

    adata = adata_utils.filter_low_quality_cells(
        adata,
        save_path=save_path / "plots",
        **({"min_counts": min_counts} if min_counts is not None else {}),
        logger=logger,
    )

    adata = adata_utils.filter_genes(
        adata, save_path=save_path / "plots", logger=logger
    )

    adata = adata_utils.normalize_counts(
        adata, save_path=save_path / "plots", seg_method=args.seg_method, logger=logger
    )

    adata = adata_utils.pca_umap_single(
        adata, save_path=save_path / "plots", logger=logger
    )

    adata = adata_utils.integration_harmony(
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


if __name__ == "__main__":
    main()
