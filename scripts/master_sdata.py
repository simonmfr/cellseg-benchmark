import logging
import os
import sys
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

sample = sys.argv[1]  # sample name
data_path = sys.argv[2]  # merscope data
zmode = sys.argv[3]
data_dir = sys.argv[4]  # base directory
if len(sys.argv) > 5:
    n_ficture = int(sys.argv[5])
else:
    n_ficture = 21

logger.info("Importing images and points")
su.process_merscope(sample, data_dir, data_path, zmode=zmode)

logger.info("Importing master sdata")
sdata_path = join(data_dir, "samples", sample)
sdata_main = read_zarr(join(sdata_path, "sdata_z3.zarr"))

# only work on methods with actual data available
seg_methods = [
    method
    for method in os.listdir(join(sdata_path, "results"))
    if os.path.isdir(join(sdata_path, "results", method, "sdata.zarr"))
]
su.integrate_segmentation_data(
    sdata_path, seg_methods, sdata_main, write_to_disk=True, data_path=data_path, logger=logger
)
