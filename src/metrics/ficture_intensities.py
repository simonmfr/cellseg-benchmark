import sys
from os import listdir
from os.path import join
from pathlib import Path

import numpy as np
import pandas as pd
from sopa.aggregation import aggregate_channels
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from tifffile import imread

sys.path.insert(1, str(Path(__file__).parent.parent.resolve()))
import ficture_utils


def ficture_intensities(
    sdata: SpatialData, data_path, ficture_path, shapes_key: str
) -> pd.DataFrame:
    """Compute average intensities for ficture factors in cells given by sdata.

    Args:
        sdata: sdata from segmentation containing at least cell boundaries.
        data_path: path to the original merscope output
        ficture_path: path to the ficture output
        shapes_key: key for valid cell boundaries

    Returns: average intensities for ficture factors in cells given by sdata.

    """
    DAPI_shape = imread(join(data_path, "images/mosaic_DAPI_z3.tif")).shape
    transform = pd.read_csv(
        join(data_path, "images/micron_to_mosaic_pixel_transform.csv"), sep=" ", header=None
    )

    ficture_full_path = ""
    for file in listdir(ficture_path):
        if file.endswith(".pixel.sorted.tsv.gz"):
            ficture_full_path = join(ficture_path, file)
    assert ficture_full_path != "", (
        "Ficture path not correct or Ficture output not computed."
    )

    fic_header = ["BLOCK", "X", "Y", "K1", "K2", "K3", "P1", "P2", "P3"]
    ficture_pixels = pd.read_csv(
        ficture_full_path, sep="\t", names=fic_header, comment="#"
    )

    metadata = ficture_utils.parse_metadata(ficture_full_path)
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

    for factor in unique_factors:
        try:
            image_stack
        except NameError:
            image_stack = ficture_utils.create_factor_level_image(ficture_pixels, factor, DAPI_shape)
        else:
            image_stack = np.concatenate(
                (
                    image_stack,
                    ficture_utils.create_factor_level_image(ficture_pixels, factor, DAPI_shape),
                ),
                axis=0, dtype=np.uint16,
            )
    sdata["image"] = Image2DModel.parse(image_stack)

    intensities = aggregate_channels(sdata, image_key="image", shapes_key=shapes_key)
    pd_intensity = pd.DataFrame(
        intensities,
        index=sdata[shapes_key].index,
        columns=[f"ficture_{i}_intensity" for i in unique_factors], #TODO: also get ficture factors (num)
    )
    pd_intensity = pd_intensity/pd_intensity.max(axis=None)
    return pd_intensity
