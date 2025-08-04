import argparse
import logging
import os
import warnings
from os.path import join

from spatialdata import read_zarr

from cellseg_benchmark import sdata_utils as su

warnings.filterwarnings("ignore")

logger = logging.getLogger("shape_mapping")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

parser = argparse.ArgumentParser(
    description="Creates the master sdata, collecting all segmentations from subfolders."
)
parser.add_argument("sample", help="sample name.")
parser.add_argument("data_path", help="Path to data folder with merscope data.")
parser.add_argument(
    "zmode", choices=["z3", "3d"], help="mode of master sdata. Either 'z3' or '3d'."
)
parser.add_argument("data_dir", help="output data folder.")
parser.add_argument(
    "--n_ficture",
    default=21,
    type=int,
    help="Consider Ficture model with n_ficture factors.",
)
args = parser.parse_args()

logger.info("Importing images and points")
su.process_merscope(args.sample, args.data_dir, args.data_path, zmode=args.zmode)

logger.info("Importing master sdata")
sdata_path = join(args.data_dir, "samples", args.sample)
sdata_main = read_zarr(join(sdata_path, "sdata_z3.zarr"))

# only work on methods with actual data available
seg_methods = [
    method
    for method in os.listdir(join(sdata_path, "results"))
    if os.path.isdir(join(sdata_path, "results", method, "sdata.zarr"))
]
su.integrate_segmentation_data(
    sdata_path,
    seg_methods,
    sdata_main,
    write_to_disk=True,
    data_path=args.data_path,
    logger=logger,
)
logger.info("Finished")
