import sys
from os.path import join
from pathlib import Path
from subprocess import run

from pandas import read_csv
from sopa.io.explorer import write
from spatialdata import read_zarr
from spatialdata_io import merscope

data_path = sys.argv[1]
save_path = sys.argv[2]

sdata = merscope(
    data_path,
    vpt_outputs={
        "cell_by_gene": Path(join(data_path, "cell_by_gene.csv")),
        "cell_metadata": Path(join(data_path, "cell_metadata.csv")),
        "cell_boundaries": Path(join(data_path, "cell_boundaries.parquet")),
    },
)
translation = read_csv(
    join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
    sep=" ",
    header=None,
)
sdata.write(join(save_path, "sdata.zarr"))
sdata = read_zarr(join(save_path, "sdata.zarr"))

write(
    join(save_path, "sdata.explorer"),
    sdata,
    gene_column="gene",
    ram_threshold_gb=4,
    pixel_size=1 / translation.iloc[0, 0],
)
run(["rm", "-r", join(save_path, "sdata.zarr", "images")])
run(["rm", "-r", join(save_path, "sdata.zarr", "points")])
