import argparse
from os.path import join
from pathlib import Path
from subprocess import run

from pandas import read_csv
from sopa.io.explorer import write
from spatialdata import read_zarr
from spatialdata_io import merscope

parser = argparse.ArgumentParser(
    description="Compute default segmentation from merscope."
)
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument("save_path", help="Path to output folder.")
parser.add_argument(
    "--explorer", action="store_true", help="if to compute explorer files"
)
args = parser.parse_args()

sdata = merscope(
    args.data_path,
    vpt_outputs={
        "cell_by_gene": Path(join(args.data_path, "cell_by_gene.csv")),
        "cell_metadata": Path(join(args.data_path, "cell_metadata.csv")),
        "cell_boundaries": Path(join(args.data_path, "cell_boundaries.parquet")),
    },
)
sdata["table"].obs.rename(columns={"EntityID": "cell_id"}, inplace=True)

boundaries = sdata[list(sdata.shapes.keys())[0]]
boundaries.rename(columns={"EntityID": "cell_id"}, inplace=True)
boundaries.index = boundaries["cell_id"]

sdata["table"].uns["spatialdata_attrs"]["instance_key"] = "cell_id"

translation = read_csv(
    join(args.data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
    sep=" ",
    header=None,
)
sdata.write(join(args.save_path, "sdata.zarr"))
sdata = read_zarr(join(args.save_path, "sdata.zarr"))

if args.explorer:
    write(
        join(args.save_path, "sdata.explorer"),
        sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )
run(["rm", "-r", join(args.save_path, "sdata.zarr", "images")])
run(["rm", "-r", join(args.save_path, "sdata.zarr", "points")])
