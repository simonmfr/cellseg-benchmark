import os
import sys
from os import listdir
from os.path import join
from pathlib import Path
from re import split

import numpy as np
import pandas as pd
from spatialdata import read_zarr
from tqdm import tqdm
from tifffile import imread

sys.path.insert(1, join(str(Path(__file__).parent.parent.resolve()), "src"))
from ficture_utils import create_factor_level_image, parse_metadata
from metrics.ficture_intensities import ficture_intensities

sys.path.insert(1, join(Path(__file__).parent.parent.resolve(), "src"))
import sdata_utils as su

sample = sys.argv[1]  # sample name
data_path = sys.argv[2]  # merscope data
zmode = sys.argv[3]
data_dir = sys.argv[4]  # base directory
if len(sys.argv) > 5:
    n_ficture = int(sys.argv[5])
else:
    n_fiture = 21

su.process_merscope(sample, data_dir, data_path, zmode=zmode)

sdata_path = join(data_dir, "samples", sample)
sdata_main = read_zarr(join(sdata_path, "sdata_z3.zarr"))

# only work on methods with actual data available
seg_methods = [
    method
    for method in os.listdir(join(sdata_path, "results"))
    if os.path.isdir(join(sdata_path, "results", method, "sdata.zarr"))
]
su.integrate_segmentation_data(
    sdata_path, seg_methods, sdata_main, write_to_disk=True, data_path=data_path
)

DAPI_shape = imread(join(data_path, "images/mosaic_DAPI_z3.tif")).shape
transform = pd.read_csv(
    join(data_path, "images/micron_to_mosaic_pixel_transform.csv"), sep=" ", header=None
)

assert "Ficture" in os.listdir(join(sdata_path, "results")), f"No Ficture computed for {sample}."
ficture_path = join(sdata_path, "results", "Ficture", "output")
for file in listdir(ficture_path):
    if n_ficture == int(split(r"\.|F", file)[1]):
        ficture_path = join(ficture_path, file)

ficture_full_path = ""
for file in listdir(ficture_path):
    if file.endswith(".pixel.sorted.tsv.gz"):
        ficture_full_path = join(ficture_path, file)
        n_factors = int(split(r"\.|F", file)[1])
assert ficture_full_path != "", (
    "Ficture output not correctly computed."
)

fic_header = ["BLOCK", "X", "Y", "K1", "K2", "K3", "P1", "P2", "P3"]
ficture_pixels = pd.read_csv(ficture_full_path, sep="\t", names=fic_header, comment="#")

metadata = parse_metadata(ficture_full_path)
scale = float(metadata["SCALE"])
offset_x = float(metadata["OFFSET_X"])
offset_y = float(metadata["OFFSET_Y"])
ficture_pixels["X_pixel"] = (
    ficture_pixels["X"] / scale * transform.iloc[0, 0]
    + offset_x * transform.iloc[0, 0]
    + transform.iloc[0, 2]
)
ficture_pixels["Y_pixel"] = (
    ficture_pixels["Y"] / scale * transform.iloc[1, 1]
    + offset_y * transform.iloc[1, 1]
    + transform.iloc[1, 2]
)
del transform, metadata

unique_factors = (
    list(np.unique(ficture_pixels["K1"]))
    + list(np.unique(ficture_pixels["K2"]))
    + list(np.unique(ficture_pixels["K3"]))
)
unique_factors = list(set(unique_factors))

for factor in tqdm(unique_factors):
    try:
        image_stack
    except NameError:
        image_stack = create_factor_level_image(ficture_pixels, factor, DAPI_shape)
    else:
        image_stack = np.concatenate(
            (
                image_stack,
                create_factor_level_image(ficture_pixels, factor, DAPI_shape),
            ),
            axis=0,
            dtype=np.uint16,
        )

shape_keys = ["_".join(key.split("_")[1:]) for key in sdata_main.shapes.keys()]
table_keys = ["_".join(key.split("_")[1:]) for key in sdata_main.tables.keys()]
full_keys = list(set(shape_keys) & set(table_keys))

for key in tqdm(full_keys):
    ficture_intensities(
        sdata_main, image_stack, key, n_factors, unique_factors, update_element=True
    )
    cell_type = pd.read_csv(join(sdata_path, "results", key, "cell_type_annotation", "adata_obs_annotated.csv"))["cell_type_final"]
    sdata_main[f"adata_{key}"].obs = sdata_main[f"adata_{key}"].obs.merge(cell_type, how="left", left_index=True, right_index=True)
    sdata_main.write_element(f"adata_{key}")
