import argparse
import os
from os.path import join
from pathlib import Path
from subprocess import run

import pandas as pd
from pandas import read_csv
from sopa.aggregation import aggregate_channels
from sopa.io.explorer import write
from sopa.utils import validated_channel_names, get_spatial_image
from spatialdata import read_zarr
from spatialdata_io import merscope

parser = argparse.ArgumentParser(
    description="Transform vpt 2D pipeline output to sdata."
)
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument("save_path", help="Path to output folder.")
parser.add_argument(
    "--explorer", action="store_true", help="if to compute explorer files"
)
args = parser.parse_args()

assert any(
    [
        "cell_by_gene.csv" in file
        for file in os.listdir(join(args.save_path, "analysis_outputs"))
    ]
), "not correctly computed"
sdata = merscope(
    args.data_path,
    vpt_outputs={
        "cell_by_gene": Path(
            join(args.save_path, "analysis_outputs", "cell_by_gene.csv")
        ),
        "cell_metadata": Path(
            join(args.save_path, "analysis_outputs", "cell_metadata.csv")
        ),
        "cell_boundaries": Path(
            join(args.save_path, "analysis_outputs", "cellpose2_micron_space.parquet")
        ),
    },
)
sdata["table"].obs.rename(columns={"EntityID": "cell_id"}, inplace=True)
shapes_key = list(sdata.shapes.keys())[0]
translation = read_csv(
    join(args.data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
    sep=" ",
    header=None,
)
boundaries = sdata[shapes_key]
boundaries.rename(columns={"EntityID": "cell_id"}, inplace=True)
boundaries.set_index("cell_id", drop=False, inplace=True)
boundaries.index = boundaries.index.rename(None)

sdata["table"].uns["spatialdata_attrs"]["instance_key"] = "cell_id"
sdata["table"].obsm['intensities'] = pd.DataFrame(
    aggregate_channels(sdata, shapes_key=shapes_key),
    columns=validated_channel_names(
        get_spatial_image(sdata, list(sdata.images.keys())[0], return_key=True)[1]
    ),
    index=sdata[shapes_key].index
)

sdata.write(join(args.save_path, "sdata.zarr"), overwrite=True)
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
