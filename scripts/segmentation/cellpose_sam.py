#!/usr/bin/env python
import argparse
import pathlib
import subprocess

import pandas as pd
import sopa
import torch
from cellpose import models
from spatialdata import read_zarr
from spatialdata.transformations import get_transformation, remove_transformation

parser = argparse.ArgumentParser(description="Compute Cellpose-SAM segmentation.")
parser.add_argument("data_path", help="Path to merfish output folder.")
parser.add_argument("save_path", help="Path to output folder.")
parser.add_argument(
    "--staining",
    help="Staining used in addition to DAPI, e.g. PolyT. Default DAPI only.",
)
parser.add_argument("--model", default="cpsam_v2", help="Cellpose pretrained model.")
args = parser.parse_args()


def main(data_path, save_path, staining, model_name):
    """Cellpose-SAM (cellpose>=4) wrapped as a sopa custom staining-based segmentation."""
    assert torch.cuda.is_available(), "Cellpose-SAM requires a GPU."
    model = models.CellposeModel(gpu=True, pretrained_model=model_name)

    sdata = sopa.io.merscope(data_path)
    translation = pd.read_csv(
        pathlib.Path(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )

    sdata.write(pathlib.Path(save_path, "sdata_tmp.zarr"), overwrite=True)
    sdata = read_zarr(pathlib.Path(save_path, "sdata_tmp.zarr"))

    image_key = list(sdata.images.keys())[0]
    if "micron" in get_transformation(sdata.images[image_key], get_all=True):
        remove_transformation(sdata.images[image_key], to_coordinate_system="micron")

    sopa.make_image_patches(sdata, patch_width=8900, patch_overlap=178)
    sopa.settings.parallelization_backend = None

    sopa.segmentation.custom_staining_based(
        sdata,
        method=lambda img: model.eval(
            img, channel_axis=0, diameter=89 if staining else 60
        )[0],
        channels=[staining, "DAPI"] if staining else ["DAPI"],
        min_area=2000,
    )

    sopa.aggregate(
        sdata, gene_column="gene", aggregate_channels=True, min_transcripts=10
    )
    sopa.io.explorer.write(
        pathlib.Path(save_path, "sdata.explorer"),
        sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )

    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(pathlib.Path(save_path, "sdata.zarr"), overwrite=True)
    subprocess.run(["rm", "-r", pathlib.Path(save_path, "sdata_tmp.zarr")])


if __name__ == "__main__":
    main(args.data_path, args.save_path, args.staining, args.model)
