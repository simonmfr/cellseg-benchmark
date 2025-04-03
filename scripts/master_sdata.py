import os
import sys
from os.path import join
from pathlib import Path

from spatialdata import read_zarr

sys.path.insert(1, join(Path(__file__).parent.parent.resolve(), "src"))
import sdata_utils as su

sample = sys.argv[1]  # sample name
data_path = sys.argv[2]  # merscope data
zmode = sys.argv[3]
data_dir = sys.argv[4]  # base directory

su.process_merscope(sample, data_dir, data_path, zmode=zmode)

sdata_path = join(data_dir, "samples", sample)
sdata_main = read_zarr(join(sdata_path, "sdata_z3.zarr"))

# only work on methods with actual data available
seg_methods = [
    method
    for method in os.listdir(join(sdata_path, "results"))
    if os.path.isdir(join(sdata_path, "results", method, "sdata.zarr"))
]
su.integrate_segmentation_data(
    sdata_path, seg_methods, sdata_main, write_to_disk=True, data_path=data_path
)
