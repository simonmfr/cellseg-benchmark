#!/usr/bin/env python
import argparse
import pathlib
import subprocess

import numpy as np
import pandas as pd
import sopa
from cellSAM import cellsam_pipeline
from spatialdata import read_zarr
from spatialdata.transformations import get_transformation, remove_transformation

parser = argparse.ArgumentParser(description="Compute CellSAM segmentation.")
parser.add_argument("data_path", help="Path to merfish output folder.")
parser.add_argument("save_path", help="Path to output folder.")
parser.add_argument(
    "--use_polyt", action="store_true", help="Use DAPI + PolyT instead of DAPI only."
)
parser.add_argument("--block_size", type=int, default=512, help="CellSAM block size.")
parser.add_argument("--patch_width", type=int, default=8192, help="Sopa patch width.")
parser.add_argument(
    "--patch_overlap", type=int, default=200, help="Sopa patch overlap."
)
parser.add_argument("--min_area", type=int, default=2000, help="Minimum cell area.")
args = parser.parse_args()


def cellsam_for_sopa(image, block_size=512, use_polyt=False):
    """Run cellsam_pipeline on a sopa patch and return an integer label mask."""
    if use_polyt:
        blank = np.zeros_like(image[0:1])
        img = np.transpose(np.concatenate([blank, image], axis=0), (1, 2, 0))
    else:
        img = image[0]

    if img.dtype == np.uint16:
        img = img.astype(np.float32) / 65535.0
    elif img.dtype == np.uint8:
        img = img.astype(np.float32) / 255.0
    else:
        img = img.astype(np.float32)

    mask = cellsam_pipeline(img, use_wsi=True, block_size=block_size)
    return mask if mask is not None else np.zeros(image.shape[1:], dtype=np.int32)


def main(
    data_path, save_path, use_polyt, block_size, patch_width, patch_overlap, min_area
):
    """CellSAM algorithm wrapped as a sopa custom staining-based segmentation."""
    sdata = sopa.io.merscope(data_path)
    translation = pd.read_csv(
        pathlib.Path(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )

    sdata.write(pathlib.Path(save_path, "sdata_tmp.zarr"), overwrite=True)
    sdata = read_zarr(pathlib.Path(save_path, "sdata_tmp.zarr"))

    # remove 3D micron transform - incompatible with sopa 2D patching
    image_key = list(sdata.images.keys())[0]
    if "micron" in get_transformation(sdata.images[image_key], get_all=True):
        remove_transformation(sdata.images[image_key], to_coordinate_system="micron")

    sopa.make_image_patches(sdata, patch_width=patch_width, patch_overlap=patch_overlap)
    sopa.settings.parallelization_backend = None  # CellSAM uses the GPU exclusively

    sopa.segmentation.custom_staining_based(
        sdata,
        method=lambda img: cellsam_for_sopa(
            img, block_size=block_size, use_polyt=use_polyt
        ),
        channels=["DAPI", "PolyT"] if use_polyt else ["DAPI"],
        min_area=min_area,
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
    main(
        args.data_path,
        args.save_path,
        args.use_polyt,
        args.block_size,
        args.patch_width,
        args.patch_overlap,
        args.min_area,
    )
