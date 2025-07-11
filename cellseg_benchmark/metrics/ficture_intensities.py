from typing import Tuple, Union

import dask
import geopandas as gpd
import numpy as np
import numpy.ma as ma
import pandas as pd
import shapely
from dask.diagnostics import ProgressBar
from shapely.geometry import Polygon, box
from sopa.segmentation.shapes import expand_radius, pixel_outer_bounds, rasterize
from sopa.utils import get_boundaries, get_spatial_image, to_intrinsic
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, ShapesModel
from xarray import DataArray


def ficture_intensities(
    sdata: SpatialData,
    images: np.ndarray,
    key: str,
    n_factors: int,
    unique_factors: list,
    var: bool = True,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, pd.DataFrame]]:
    """Compute the mean and optionally the variance of the ficture intensities per cell.

    Args:
        sdata: sdata of method with "table" table
        images: ficure images in numpy stack
        key: key of method
        n_factors: number of factors in ficture model
        unique_factors: number of unique factors in top 3
        var: if variance should be calculated

    Returns: None, if update_element is True, otherwise mean and variance of ficture intensities.

    """
    boundary_key = sdata["table"].uns["spatialdata_attrs"]["region"]
    tmp_sdata = SpatialData()
    tmp_sdata["ficture_images"] = Image2DModel.parse(images)
    tmp_sdata[f"boundaries_{key}"] = ShapesModel.parse(sdata[boundary_key])

    intensities = aggregate_channels(
        tmp_sdata, image_key="ficture_images", shapes_key=f"boundaries_{key}"
    )
    pd_intensity = pd.DataFrame(
        intensities,
        index=sdata["table"].obs_names,
        columns=[f"fictureF{n_factors}_{i}_mean_intensity" for i in unique_factors],
    )
    pd_intensity = pd_intensity / pd_intensity.max(axis=None)
    if var:
        variance = aggregate_channels(
            tmp_sdata,
            image_key="ficture_images",
            shapes_key=f"boundaries_{key}",
            mode="variance",
            means=intensities,
        )
        pd_variance = pd.DataFrame(
            variance,
            index=sdata["table"].obs_names,
            columns=[f"fictureF{n_factors}_{i}_var_intensity" for i in unique_factors],
        )
        pd_variance = pd_variance / np.power(pd_intensity.max(axis=None), 2)
        return (pd_intensity, pd_variance)
    return pd_intensity


# from sopa.aggregation.channels.py
AVAILABLE_MODES = ["average", "min", "max", "variance", "sum"]


def aggregate_channels(
    sdata: SpatialData,
    image_key: str | None = None,
    shapes_key: str | None = None,
    expand_radius_ratio: float = 0,
    mode: str = "average",
    means: np.ndarray | None = None,
) -> np.ndarray:
    """Aggregate the channel intensities per cell (either `"average"`, or take the `"min"` / `"max"`).

    Args:
        sdata: A `SpatialData` object
        image_key: Key of `sdata` containing the image. If only one `images` element, this does not have to be provided.
        shapes_key: Key of `sdata` containing the cell boundaries. If only one `shapes` element, this does not have to be provided.
        expand_radius_ratio: Cells polygons will be expanded by `expand_radius_ratio * mean_radius`. This help better aggregate boundary stainings.
        mode: Aggregation mode. One of `"average"`, `"min"`, `"max"`. By default, average the intensity inside the cell mask.
        means: If given, provides the means of the channel intensities.

    Returns:
        A numpy `ndarray` of shape `(n_cells, n_channels)`
    """
    assert mode in AVAILABLE_MODES, (
        f"Invalid {mode=}. Available modes are {AVAILABLE_MODES}"
    )
    if mode == "variance":
        assert means is not None, "means required for variance computation"

    image = get_spatial_image(sdata, image_key)

    geo_df = get_boundaries(sdata, key=shapes_key)
    geo_df = to_intrinsic(sdata, geo_df, image)
    geo_df = expand_radius(geo_df, expand_radius_ratio)

    return _aggregate_channels_aligned(image, geo_df, mode, means)


def _aggregate_channels_aligned(
    image: DataArray,
    geo_df: gpd.GeoDataFrame | list[Polygon],
    mode: str,
    means: np.ndarray | None = None,
) -> np.ndarray:
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
                masked_image = ma.masked_array(
                    sub_image, 1 - np.repeat(mask[None], n_channels, axis=0)
                )
                aggregation[index] = np.minimum(
                    aggregation[index], masked_image.min(axis=(1, 2))
                )
            elif mode == "variance":
                func = np.sum
                values = func(
                    np.power(
                        sub_image * mask - means[index][:, np.newaxis, np.newaxis], 2
                    ),
                    axis=(1, 2),
                )
                aggregation[index] += values
            elif mode in ["average", "max", "sum"]:
                if mode in ["average", "sum"]:
                    func = np.sum
                else:
                    func = np.max
                values = func(sub_image * mask, axis=(1, 2))

                match mode:
                    case "average":
                        aggregation[index] += values
                    case "max":
                        aggregation[index] = np.maximum(aggregation[index], values)
                    case "sum":
                        aggregation[index] += values

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
        case "sum":
            return aggregation
