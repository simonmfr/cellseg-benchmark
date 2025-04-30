import sys
from os import listdir
from os.path import join
from pathlib import Path
from re import split
from re import search

import dask
import geopandas as gpd
import numpy as np
import numpy.ma as ma
import pandas as pd
import shapely
from dask.diagnostics import ProgressBar
from shapely.geometry import Polygon, box
from sopa.utils import get_spatial_image, get_boundaries, to_intrinsic
from sopa.segmentation.shapes import expand_radius, pixel_outer_bounds, rasterize
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from tifffile import imread
from xarray import DataArray

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
        join(data_path, "images/micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )

    ficture_full_path = ""
    for file in listdir(ficture_path):
        if file.endswith(".pixel.sorted.tsv.gz"):
            ficture_full_path = join(ficture_path, file)
            n_factors = int(search(r'nF(\d+)', file).group(1))
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
            image_stack = ficture_utils.create_factor_level_image(
                ficture_pixels, factor, DAPI_shape
            )
        else:
            image_stack = np.concatenate(
                (
                    image_stack,
                    ficture_utils.create_factor_level_image(
                        ficture_pixels, factor, DAPI_shape
                    ),
                ),
                axis=0,
                dtype=np.uint16,
            )
    sdata["image"] = Image2DModel.parse(image_stack)

    intensities = aggregate_channels(sdata, image_key="image", shapes_key=shapes_key)
    variance = aggregate_channels(sdata, image_key="image", shapes_key=shapes_key, mode="variance", means=intensities)
    pd_intensity = pd.DataFrame(
        intensities,
        index=sdata[shapes_key].index,
        columns=[f"fictureF{n_factors}_{i}_mean_intensity" for i in unique_factors],
    )
    pd_variance = pd.DataFrame(
        variance,
        index=sdata[shapes_key].index,
        columns=[f"fictureF{n_factors}_{i}_var_intensity" for i in unique_factors],
    )
    pd_intensity = pd_intensity / pd_intensity.max(axis=None)
    return pd_intensity, pd_variance

# from sopa.aggregation.channels.py
AVAILABLE_MODES = ["average", "min", "max", "variance"]

def aggregate_channels(
    sdata: SpatialData,
    image_key: str | None = None,
    shapes_key: str | None = None,
    expand_radius_ratio: float = 0,
    mode: str = "average",
    means: np.ndarray | None = None
) -> np.ndarray:
    """Aggregate the channel intensities per cell (either `"average"`, or take the `"min"` / `"max"`).

    Args:
        sdata: A `SpatialData` object
        image_key: Key of `sdata` containing the image. If only one `images` element, this does not have to be provided.
        shapes_key: Key of `sdata` containing the cell boundaries. If only one `shapes` element, this does not have to be provided.
        expand_radius_ratio: Cells polygons will be expanded by `expand_radius_ratio * mean_radius`. This help better aggregate boundary stainings.
        mode: Aggregation mode. One of `"average"`, `"min"`, `"max"`. By default, average the intensity inside the cell mask.

    Returns:
        A numpy `ndarray` of shape `(n_cells, n_channels)`
    """
    assert mode in AVAILABLE_MODES, f"Invalid {mode=}. Available modes are {AVAILABLE_MODES}"
    if mode == "variance":
        assert means is not None, "means required for variance computation"

    image = get_spatial_image(sdata, image_key)

    geo_df = get_boundaries(sdata, key=shapes_key)
    geo_df = to_intrinsic(sdata, geo_df, image)
    geo_df = expand_radius(geo_df, expand_radius_ratio)

    return _aggregate_channels_aligned(image, geo_df, mode, means)


def _aggregate_channels_aligned(image: DataArray, geo_df: gpd.GeoDataFrame | list[Polygon], mode: str, means: np.ndarray | None = None) -> np.ndarray:
    """Average channel intensities per cell. The image and cells have to be aligned, i.e. be on the same coordinate system.

    Args:
        image: A `DataArray` of shape `(n_channels, y, x)`
        geo_df: A `GeoDataFrame` whose geometries are cell boundaries (polygons)
        mode: type of aggregation
        means: A `np.ndarray` of shape `(n_channels,)` whose entries are mean intensities

    Returns:
        A numpy `ndarray` of shape `(n_cells, n_channels)`
    """
    cells = geo_df if isinstance(geo_df, list) else list(geo_df.geometry)
    tree = shapely.STRtree(cells)

    n_channels = len(image.coords["c"])
    areas = np.zeros(len(cells))
    if mode == "min":
        aggregation = np.full((len(cells), n_channels), fill_value=np.inf)
    else:
        aggregation = np.zeros((len(cells), n_channels))

    chunk_sizes = image.data.chunks
    offsets_y = np.cumsum(np.pad(chunk_sizes[1], (1, 0), "constant"))
    offsets_x = np.cumsum(np.pad(chunk_sizes[2], (1, 0), "constant"))

    def _average_chunk_inside_cells(chunk, iy, ix):
        ymin, ymax = offsets_y[iy], offsets_y[iy + 1]
        xmin, xmax = offsets_x[ix], offsets_x[ix + 1]

        patch = box(xmin, ymin, xmax, ymax)
        intersections = tree.query(patch, predicate="intersects")

        for index in intersections:
            cell = cells[index]
            bounds = pixel_outer_bounds(cell.bounds)

            sub_image = chunk[
                :,
                max(bounds[1] - ymin, 0) : bounds[3] - ymin,
                max(bounds[0] - xmin, 0) : bounds[2] - xmin,
            ]

            if sub_image.shape[1] == 0 or sub_image.shape[2] == 0:
                continue

            mask = rasterize(cell, sub_image.shape[1:], bounds)

            areas[index] += np.sum(mask)

            if mode == "min":
                masked_image = ma.masked_array(sub_image, 1 - np.repeat(mask[None], n_channels, axis=0))
                aggregation[index] = np.minimum(aggregation[index], masked_image.min(axis=(1, 2)))
            elif mode in ["average", "max", "variance"]:
                func = np.sum if mode == "average" else np.max
                values = func(sub_image * mask, axis=(1, 2))

                match mode:
                    case "average":
                        aggregation[index] += values
                    case "variance":
                        aggregation[index] += np.power(values - means[index], 2)
                    case "max":
                        aggregation[index] = np.maximum(aggregation[index], values)

    with ProgressBar():
        tasks = [
            dask.delayed(_average_chunk_inside_cells)(chunk, iy, ix)
            for iy, row in enumerate(image.chunk({"c": -1}).data.to_delayed()[0])
            for ix, chunk in enumerate(row)
        ]
        dask.compute(tasks)

    match mode:
        case "average":
            return aggregation / areas[:, None].clip(1)
        case "max":
            return aggregation
        case "variance":
            return aggregation / (areas[:, None].clip(2) - 1)
