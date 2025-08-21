import argparse
import logging
import os

import numpy as np
from pandas import read_csv
from skimage.filters import gaussian
from tifffile import imread, imwrite

parser = argparse.ArgumentParser(
    description="Compute tiffs based on transcript locations. Uses 10-times gaußian kernel with sigma=3."
)
parser.add_argument("path", help="Path to merscope output folder.")
args = parser.parse_args()

logger = logging.getLogger("annotation")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

image_dir = os.path.join(args.path, "images")
logger.info("Determine need for Transcript tiffs")
dapi_levels = {
    f.split(".")[0].split("_")[-1]
    for f in os.listdir(image_dir)
    if "DAPI" in f
}
transcript_levels = {
    f.split(".")[0].split("_")[-1]
    for f in os.listdir(image_dir)
    if "Transcripts" in f
}
z_levels = dapi_levels - transcript_levels
logger.info(f"Reading transcript files from {args.path}")
image_DAPI = imread(os.path.join(image_dir, "mosaic_DAPI_z3.tif"), key=0)
mean = np.mean(image_DAPI)
bins_y = np.linspace(0, image_DAPI.shape[1], num=image_DAPI.shape[1] + 1)
bins_x = np.linspace(0, image_DAPI.shape[0], num=image_DAPI.shape[0] + 1)
del image_DAPI

transcripts = read_csv(os.path.join(args.path, "detected_transcripts.csv"))[
    ["global_x", "global_y", "global_z", "gene"]
]
transform = read_csv(
    os.path.join(args.path, "images/micron_to_mosaic_pixel_transform.csv"),
    names=["x", "y", "z"],
    delimiter=" ",
)
transcripts = transcripts[transcripts["gene"] != "Mbp"]
transcripts["global_x"] = (
    transcripts["global_x"].values * transform["x"].values[0] + transform["z"].values[0]
)
transcripts["global_y"] = (
    transcripts["global_y"].values * transform["y"].values[1] + transform["z"].values[1]
)

del transform

for z in z_levels:
    logger.info(f"Processing z={z[1:]}")
    transcripts_slice = transcripts[transcripts["global_z"] == int(z[1:])]

    image, _, _ = np.histogram2d(
        transcripts_slice["global_y"],
        transcripts_slice["global_x"],
        bins=[bins_x, bins_y],
    )

    for _ in range(10):
        image = gaussian(image, sigma=3)

    image = (mean / np.max(image) * image).astype("uint16")
    imwrite(os.path.join(image_dir, f"mosaic_Transcripts_{z}.tif"), image)
logger.info(f"Done")