
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
from spatialdata.models import Image2DModel
from xarray import DataArray


def ficture_intensities(
    master_sdata: SpatialData,
    images,
    key,
    n_factors,
    unique_factors,
    var: bool = True,
    update_element: bool = False,
):
    """Compute the mean and optionally the variance of the ficture intensities per cell.

    Args:
        master_sdata: master sdata created by master_sdata.py
        images: ficure images in numpy stack
        key: key of method
        n_factors: number of factors in ficture model
        unique_factors: number of unique factors in top 3
        var: if variance should be calculated
        update_element: if the master sdata should be updated or not

    Returns: None, if update_element is True, otherwise mean and variance of ficture intensities.

    """
    master_sdata["ficture_images"] = Image2DModel.parse(images)
    intensities = aggregate_channels(
        master_sdata, image_key="ficture_images", shapes_key=f"boundaries_{key}"
    )
    pd_intensity = pd.DataFrame(
        intensities,
        index=master_sdata[f"adata_{key}"].obs.index,
        columns=[f"fictureF{n_factors}_{i}_mean_intensity" for i in unique_factors],
    )
    pd_intensity = pd_intensity / pd_intensity.max(axis=None)
    if var:
        variance = aggregate_channels(
            master_sdata,
            image_key="ficture_images",
            shapes_key=f"boundaries_{key}",
            mode="variance",
            means=intensities,
        )
        pd_variance = pd.DataFrame(
            variance,
            index=master_sdata[f"adata_{key}"].obs.index,
            columns=[f"fictureF{n_factors}_{i}_var_intensity" for i in unique_factors],
        )
    else:
        pd_variance = None
    if update_element:
        master_sdata[f"adata_{key}"].obsm["mean_intensity"] = pd_intensity
        if var:
            master_sdata[f"adata_{key}"].obsm["variance_intensity"] = pd_variance
        master_sdata.write_element(f"adata_{key}")
        return
    return pd_intensity, pd_variance


# from sopa.aggregation.channels.py
AVAILABLE_MODES = ["average", "min", "max", "variance"]


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
