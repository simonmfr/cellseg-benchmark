import logging
import os

import dask
import dask.array as da
import dask.diagnostics
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1
import numpy as np
import numpy.ma as ma
import ovrlpy
import pandas as pd
import polars
import shapely
import sopa.segmentation.shapes
import xarray


def compute_ovrl(
    sample: str, sample_dir: str, data_dir: str, logger: logging.Logger = None, **kwargs
) -> None:
    """Compute ovrlpy output.

    Args:
        sample: name of sample
        sample_dir: Sample directory
        data_dir: Directory to transcript .csv.
        logger: logging.Logger instance.

    Returns:
        None
    """
    if not os.path.exists(
        os.path.join(sample_dir, "vertical_doublets_ovrlpy_output.npz")
    ):
        coords_df = pd.read_csv(
            os.path.join(data_dir, "detected_transcripts.csv"), index_col=0
        )[["gene", "x", "y", "global_z"]].rename(columns={"global_z": "z"})
        run_ovrlpy(sample, coords_df, sample_dir)
    else:
        if logger:
            logger.warning(f"Ovrlpy output exists, skipping ovrlpy run for {sample}")
        else:
            print("Ovrlpy output exists, skipping ovrlpy run.")


def run_ovrlpy(
    sample_name: str,
    coordinate_df: pd.DataFrame | polars.DataFrame,
    data_dir: str,
    logger: logging.Logger = None,
) -> None:
    """Processes MERSCOPE sample using ovrlpy to detect vertical doublets.

    Reads transcript coordinates, runs ovrlpy analysis,
    and saves the integrity and signal matrices as a compressed .npz file.

    Args:
        sample_name: Name of the MERSCOPE sample.
        coordinate_df: Coordinate DataFrame of transcript coordinates.
        data_dir: Base directory to save output files.
        logger: logging.Logger instance.

    Returns:
        None
    """
    n_workers = os.cpu_count()

    if logger:
        logger.info(f"[{sample_name}] Starting analysis...")
    else:
        print(f"[{sample_name}] Starting analysis...")

    ovrlpy_obj = ovrlpy.Ovrlp(coordinate_df, n_workers=n_workers, random_state=42)
    ovrlpy_obj.analyse()

    output_path = os.path.join(data_dir, "vertical_doublets_ovrlpy_output.npz")
    np.savez_compressed(
        output_path, integrity=ovrlpy_obj.integrity_map, signal=ovrlpy_obj.signal_map
    )

    if logger:
        logger.info(
            f"[{sample_name}] Analysis complete. Results saved to: data_dir/vertical_doublets_ovrlpy_output.npz"
        )
    else:
        print(
            f"[{sample_name}] Analysis complete. Results saved to: data_dir/vertical_doublets_ovrlpy_output.npz"
        )


def compute_mean_vsi_per_polygon(
    integrity_map: np.ndarray,
    boundaries: gpd.GeoDataFrame,
    transform_matrix: np.ndarray,
    **kwargs,
) -> pd.DataFrame:
    """Compute mean vsi per polygon.

    Args:
        integrity_map: integrity map from ovrlpy. Assumed 2D
        boundaries: boundaries in microns
        transform_matrix: transformation matrix of MERSCOPE

    Returns:
        Mean vsi per polygon.
    """

    def micron_to_pixel_coords(
        geom: shapely.geometry.base.BaseGeometry,
        transform_matrix: np.ndarray,
        pixel_offset=(13, 13),
    ) -> shapely.geometry.base.BaseGeometry:
        """Transform micron coordinates to ovrlpy coordinates."""
        sx, sy = transform_matrix[0, 0], transform_matrix[1, 1]
        ox, oy = transform_matrix[0, 2], transform_matrix[1, 2]
        omx, omy = ox / sx, oy / sy

        geom = shapely.affinity.affine_transform(geom, [1, 0, 0, 1, omx, omy])
        geom = shapely.affinity.translate(
            geom, xoff=-pixel_offset[0], yoff=-pixel_offset[1]
        )
        return geom

    pic = np.expand_dims(integrity_map, axis=0)
    pic = da.from_array(pic)
    pic = xarray.DataArray(pic, dims=["c", "y", "x"])

    boundaries["geometry"] = boundaries.geometry.apply(
        micron_to_pixel_coords, args=(transform_matrix, (13, 13))
    )

    result = aggregate_channels_aligned(pic, boundaries, "average")
    var = aggregate_channels_aligned(pic, boundaries, "variance", means=result)
    result = pd.DataFrame(result, index=boundaries.index, columns=["mean_integrity"])
    result["variance"] = var
    return result


def plot_vsi_overview(
    integrity_map,
    signal_map,
    boundaries_aligned,
    vsi_mean,
    sample_name,
    boxes=None,
    png_path=None,
):
    """Plot vsi per polygon.

    Args:
        integrity_map: ovrlpy integrity map
        signal_map: ovrlpy signal map
        boundaries_aligned: boundaries in micron coordinates
        vsi_mean: mean vsi per polygon
        sample_name: sample name
        boxes: Box boundaries
        png_path: save path for image

    Returns:
        None, saves figure if png_path provided.
    """
    ny, nx = integrity_map.shape
    if boxes is None:
        boxes = [
            (int(2 / 8 * nx), int(3 / 8 * ny), int(1 / 8 * nx), int(1 / 8 * ny)),
            (int(6 / 8 * nx), int(5 / 8 * ny), int(1 / 8 * nx), int(1 / 8 * ny)),
        ]

    plot_kwargs = {"cmap": "coolwarm_r", "origin": "lower", "vmin": 0, "vmax": 1}
    boundary_kwargs = {"facecolor": "none", "edgecolor": "black"}
    _SIGNAL_THRESHOLD = 2  # from ovrlpy source code
    alpha = (signal_map / _SIGNAL_THRESHOLD).clip(
        0, 1
    ) ** 2  # fade out for pixels with low transcript signal

    fig, axs = plt.subplots(2, 2, figsize=(12, 10))

    # overview
    im = axs[0, 0].imshow(integrity_map, alpha=alpha, **plot_kwargs)
    boundaries_aligned.plot(
        ax=axs[0, 0], aspect="equal", linewidth=0.1, **boundary_kwargs
    )
    axs[0, 0].set_title(sample_name)

    for x, y, w, h in boxes:
        rect = mpl.patches.Rectangle(
            (x, y), w, h, linewidth=1, edgecolor="red", facecolor="none", linestyle="--"
        )
        axs[0, 0].add_patch(rect)

    cax = mpl_toolkits.axes_grid1.make_axes_locatable(axs[0, 0]).append_axes(
        "right", size="4%", pad=0.05
    )
    fig.colorbar(im, cax=cax).ax.set_title("VSI", fontsize=10)

    # histogram
    axs[0, 1].hist(vsi_mean, bins=100, zorder=3, color="slategrey")
    axs[0, 1].axvline(0.5, color="gray", linestyle="--", linewidth=1)
    axs[0, 1].grid(True, zorder=0)
    axs[0, 1].set_title("Mean VSI per cell - " + sample_name)

    # detailed views
    for i, (x, y, w, h) in enumerate(boxes):
        ax = axs[1, i]
        ax.imshow(integrity_map, alpha=alpha, **plot_kwargs)
        boundaries_aligned.plot(ax=ax, linewidth=0.5, **boundary_kwargs)
        ax.set_xlim(x, x + w)
        ax.set_ylim(y, y + h)
        ax.set_title(f"Box {chr(65 + i)}")

    plt.tight_layout()

    if png_path:
        plt.savefig(png_path, dpi=200)
    plt.close()


# from sopa.aggregation.channels.py
AVAILABLE_MODES = ["average", "min", "max", "variance", "sum"]


def aggregate_channels_aligned(
    image: xarray.DataArray,
    geo_df: gpd.GeoDataFrame | list[shapely.geometry.Polygon],
    mode: str,
    means: np.ndarray | None = None,
) -> np.ndarray:
    """Reduce each image channel to one value per cell polygon.

    ``mode`` selects the per-cell statistic: average, min, max, variance, or sum.
    Image and polygons must share a coordinate system. Generic image aggregation
    (e.g. ovrlpy intensities), not FICTURE-specific.

    Args:
        image: Image ``DataArray`` of shape ``(n_channels, y, x)``.
        geo_df: Cell boundary polygons (``GeoDataFrame`` or list of polygons).
        mode: Aggregation statistic: "average", "min", "max", "variance", or "sum".
        means: Per-channel means, required only when ``mode="variance"``.

    Returns:
        Array of shape ``(n_cells, n_channels)``.
    """
    assert mode in AVAILABLE_MODES, (
        f"Invalid {mode=}. Available modes are {AVAILABLE_MODES}"
    )
    if mode == "variance":
        assert means is not None, "means required for variance computation"

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

        patch = shapely.geometry.box(xmin, ymin, xmax, ymax)
        intersections = tree.query(patch, predicate="intersects")

        for index in intersections:
            cell = cells[index]
            bounds = sopa.segmentation.shapes.pixel_outer_bounds(cell.bounds)

            sub_image = chunk[
                :,
                max(bounds[1] - ymin, 0) : bounds[3] - ymin,
                max(bounds[0] - xmin, 0) : bounds[2] - xmin,
            ]

            if sub_image.shape[1] == 0 or sub_image.shape[2] == 0:
                continue

            mask = sopa.segmentation.shapes.rasterize(cell, sub_image.shape[1:], bounds)

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

    with dask.diagnostics.ProgressBar():
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
