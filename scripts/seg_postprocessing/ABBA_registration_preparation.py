import argparse
import logging
import warnings
from pathlib import Path

import cellseg_benchmark
import pandas as pd
import spatialdata

warnings.filterwarnings("ignore")

logger = logging.getLogger("ABBA_registration_preparation")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

def main():
    parser = argparse.ArgumentParser(
        description=("Create points table for registration in QuPath.")
    )
    parser.add_argument("sample", type=str, help="Sample name.")
    parser.add_argument("out_dir", help="Path to output directory.")
    parser.add_argument("base_path", default=cellseg_benchmark.BASE_PATH, help="base path of data structure.")
    args = parser.parse_args()

    zarr_path = Path(args.base_path) / "samples" / args.sample / "sdata_z3.zarr"
    logger.info(f"Read zarr from {zarr_path}")
    sdata = spatialdata.read_zarr(zarr_path, selection=("tables",))

    all_points = []
    for key, value in sdata.tables.items():
        points = pd.DataFrame(value.obsm['spatial'], index=value.obs_names.values, columns=["x", "y"])
        points['cell_id'] = value.obs_names
        points['sample'] = args.sample
        points['segmentation'] = "_".join(key.split("_")[1:])
        points = points.loc[:,['cell_id', "x", "y", "sample", "segmentation"]]
        all_points.append(points)
    all_points = pd.concat(all_points)
    save_path = Path(args.out_dir) / f"{args.sample}_points.csv"
    logger.info(f"Saving points to {save_path}")
    all_points.to_csv(save_path)