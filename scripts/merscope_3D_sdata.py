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
    transcripts=False,
    mosaic_images=False,
    vpt_outputs={
        "cell_by_gene": Path(join(save_path, "analysis_outputs", "cell_by_gene.csv")),
        "cell_metadata": Path(join(save_path, "analysis_outputs", "cell_metadata.csv")),
        "cell_boundaries": Path(join(save_path, "analysis_outputs", "cellpose2_micron_space.parquet")),
    },
    z_layers=[0,1,2,3,4,5,6]
)

sdata.write(join(save_path, "sdata.zarr"))
