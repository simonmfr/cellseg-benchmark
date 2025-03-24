import os
import sys

import numpy as np
from pandas import read_csv
from skimage.filters import gaussian
from tifffile import imread, imwrite

path = sys.argv[1]

image_DAPI = imread(os.path.join(path, "images/mosaic_DAPI_z3.tif"), key=0)
mean = np.mean(image_DAPI)
bins_y = np.linspace(0, image_DAPI.shape[1], num=image_DAPI.shape[1] + 1)
bins_x = np.linspace(0, image_DAPI.shape[0], num=image_DAPI.shape[0] + 1)
del image_DAPI

transcripts = read_csv(os.path.join(path, "detected_transcripts.csv"))[
    ["global_x", "global_y", "global_z", "gene"]
]
transform = read_csv(
    os.path.join(path, "images/micron_to_mosaic_pixel_transform.csv"),
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

for i in range(7):
    transcripts_slice = transcripts[transcripts["global_z"] == i]

    image, _, _ = np.histogram2d(
        transcripts_slice["global_y"],
        transcripts_slice["global_x"],
        bins=[bins_x, bins_y],
    )

    for _ in range(10):
        image = gaussian(image, sigma=3)

    image = (mean / np.max(image) * image).astype("uint16")
    imwrite(os.path.join(path, f"images/mosaic_Transcripts_z{i}.tif"), image)
