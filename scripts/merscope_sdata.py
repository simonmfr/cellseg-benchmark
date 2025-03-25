import sys
from os.path import join
from subprocess import run

from sopa.io.explorer import write
from spatialdata_io import merscope

data_path = sys.argv[1]
save_path = sys.argv[2]

sdata = merscope(data_path)
write(join(save_path, "sdata.explorer", sdata))
sdata.write(join(save_path, "sdata.zarr"))
run(["rm", "-r", join(save_path, "sdata.zarr", "images")])
run(["rm", "-r", join(save_path, "sdata.zarr", "points")])
