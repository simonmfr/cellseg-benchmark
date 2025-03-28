import sys
from os.path import join
from pathlib import Path
from subprocess import run

from sopa.io.explorer import write
from spatialdata_io import merscope
from spatialdata import read_zarr

data_path = sys.argv[1]
save_path = sys.argv[2]

sdata = merscope(data_path,
                 vpt_outputs={'cell_by_gene': Path(join(data_path, 'cell_by_gene.csv')),
                              'cell_metadata': Path(join(data_path, 'cell_metadata.csv')),
                              'cell_boundaries': Path(join(data_path, 'cell_boundaries.parquet'))}
                 )
sdata.write(join(save_path, "sdata.zarr"))
sdata = read_zarr(join(save_path, "sdata.zarr"))

write(join(save_path, "sdata.explorer"), sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=0.108)
run(["rm", "-r", join(save_path, "sdata.zarr", "images")])
run(["rm", "-r", join(save_path, "sdata.zarr", "points")])
