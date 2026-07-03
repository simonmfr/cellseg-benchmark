#!/usr/bin/env python
import argparse
import logging
import pathlib

import spatialdata as sd

import cellseg_benchmark.sdata_utils as sdata_utils

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s]: %(message)s"
)
logger = logging.getLogger("shape_mapping")


def main():
    p = argparse.ArgumentParser(
        description="Creates a master sdata for a given sample, containing multiple segmentation results."
    )
    p.add_argument("sample", help="Sample name.")
    p.add_argument(
        "data_path",
        help="Path to folder with merscope output data (e.g. /cohort1/slide2/region0).",
    )
    p.add_argument(
        "zmode",
        choices=["z3"],
        help="Mode of master sdata. Either 'z3' or '3d' (currently only z3 is implemented).",
    )
    p.add_argument("data_dir", help="Output data folder.")
    p.add_argument(
        "--n_ficture",
        default=21,
        type=int,
        help="Consider Ficture model with n_ficture factors.",
    )
    p.add_argument("--run_date", type=str, help="run date (YYYYMMDD).", default=None)
    p.add_argument("--organism", type=str, help="organism.", default=None)
    p.add_argument("--slide", type=str, help="slide.", default=None)
    p.add_argument("--region", type=str, help="region.", default=None)
    p.add_argument("--cohort", type=str, help="cohort.", default=None)
    p.add_argument(
        "--obs",
        action="append",
        default=[],
        metavar="KEY=VAL",
        help="Extra covariates to add to adata.obs (repeatable), e.g. --obs tissue=brain.",
    )
    args = p.parse_args()

    NONES = {"", "None", "none", "null", "NULL", None}
    for k in ["organism", "slide", "region", "cohort"]:
        if getattr(args, k) in NONES:
            setattr(args, k, None)

    extra_obs = {}
    for kv in args.obs:
        k, v = kv.split("=", 1)
        extra_obs[k] = None if v in NONES else v

    logger.info("Importing images and points...")
    sdata_utils.process_merscope(
        args.sample, args.data_dir, args.data_path, zmode=args.zmode
    )

    sdata_path = pathlib.Path(args.data_dir) / "samples" / args.sample
    sdata_main = sd.read_zarr(str(sdata_path / "sdata_z3.zarr"))

    logger.info(
        "Integrating segmentation data from available methods into main sdata..."
    )
    results_path = sdata_path / "results"
    seg_methods = [
        method.name
        for method in results_path.iterdir()
        if (method / "sdata.zarr").is_dir()
    ]
    sdata_utils.integrate_segmentation_data(
        str(sdata_path),
        seg_methods,
        sdata_main,
        run_date=args.run_date,
        organism=args.organism,
        slide=args.slide,
        region=args.region,
        cohort=args.cohort,
        write_to_disk=True,
        data_path=args.data_path,
        logger=logger,
        **extra_obs,
    )
    logger.info("Done.")


if __name__ == "__main__":
    main()
