import logging
import os
from os.path import exists, join
from typing import List, Tuple

import pandas as pd
from dask.array import from_array
from geopandas import GeoDataFrame
from matplotlib import patches, pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely import affinity
from shapely.geometry.base import BaseGeometry
from spatialdata import SpatialData
from tqdm import tqdm
import ovrlpy
import numpy as np
from polars import DataFrame
from xarray import DataArray

from .ficture_intensities import _aggregate_channels_aligned


def compute_ovrl(sdata_list: List[Tuple[str, SpatialData]], data_dir: str, logger: logging.Logger=None) -> None:
    """
    Compute ovrlpy output.
    Args:
        sdata_list: List of sdatas with names
        data_dir: Base directory for output files.
        logger: logging.Logger instance.

    Returns:
        None
    """
    for name, sdata in tqdm(sdata_list):
        if not exists(join(data_dir, "samples", name, "vertical_doublets_ovrlpy_output.npz")):
            coords_df = sdata.points[name + "_transcripts"][["gene", "x", "y", "global_z"]].compute().rename(
                columns={"global_z": "z"})
            run_ovrlpy(name, coords_df, data_dir)
        else:
            if logger:
                logger.warning(f"Ovrlpy output exists, skipping ovrlpy run for {name}")
            else:
                print(f"Ovrlpy output exists, skipping ovrlpy run.")


def run_ovrlpy(sample_name: str, coordinate_df: pd.DataFrame | DataFrame, data_dir: str, logger: logging.Logger=None) -> None:
    """
    Processes MERSCOPE sample using ovrlpy to detect vertical doublets.

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

    output_path = os.path.join(data_dir, "samples", sample_name, "vertical_doublets_ovrlpy_output.npz")
    np.savez_compressed(output_path,
                        integrity=ovrlpy_obj.integrity_map,
                        signal=ovrlpy_obj.signal_map
                        )

    if logger:
        logger.info(f"[{sample_name}] Analysis complete. Results saved to: data_dir/{os.path.join('samples', sample_name, 'vertical_doublets_ovrlpy_output.npz')}")
    else:
        print(
        f"[{sample_name}] Analysis complete. Results saved to: data_dir/{os.path.join('samples', sample_name, 'vertical_doublets_ovrlpy_output.npz')}")

def compute_mean_vsi_per_polygon(integrity_map: np.ndarray, boundaries: GeoDataFrame, transform_matrix: np.ndarray) -> pd.DataFrame:
    """
    Compute mean vsi per polygon.
    Args:
        integrity_map: integrity map from ovrlpy. Assumed 2D
        boundaries: boundaries in microns
        transform_matrix: transformation matrix of MERSCOPE

    Returns:
        Mean vsi per polygon.
    """
    def micron_to_pixel_coords(
            geom: BaseGeometry,
            transform_matrix: np.ndarray,
            pixel_offset=(13, 13),
    ) -> BaseGeometry:
        """Transform micron coordinates to ovrlpy coordinates."""
        sx, sy = transform_matrix[0, 0], transform_matrix[1, 1]
        ox, oy = transform_matrix[0, 2], transform_matrix[1, 2]
        omx, omy = ox / sx, oy / sy

        geom = affinity.affine_transform(geom, [1, 0, 0, 1, omx, omy])
        geom = affinity.translate(geom, xoff=-pixel_offset[0], yoff=-pixel_offset[1])
        return geom

    pic = np.expand_dims(integrity_map, axis=0)
    pic = from_array(pic)
    pic = DataArray(pic, dims=["c", "y", "x"])

    boundaries["geometry"] = boundaries.geometry.apply(micron_to_pixel_coords, args=(transform_matrix, (13,13)))

    result = _aggregate_channels_aligned(pic, boundaries, "average")
    return pd.DataFrame(result, index=boundaries.index, columns=["mean_integrity"])


def plot_vsi_overview(integrity_map, signal_map, boundaries_aligned, vsi_mean, sample_name, boxes=None, png_path=None):
    ny, nx = integrity_map.shape
    if boxes is None:
        boxes = [(int(2 / 8 * nx), int(3 / 8 * ny), int(1 / 8 * nx), int(1 / 8 * ny)),
                 (int(6 / 8 * nx), int(5 / 8 * ny), int(1 / 8 * nx), int(1 / 8 * ny))]

    plot_kwargs = {'cmap': 'coolwarm_r', 'origin': 'lower', 'vmin': 0, 'vmax': 1}
    boundary_kwargs = {'facecolor': 'none', 'edgecolor': 'black'}
    _SIGNAL_THRESHOLD = 2  # from ovrlpy source code
    alpha = (signal_map / _SIGNAL_THRESHOLD).clip(0, 1) ** 2  # fade out for pixels with low transcript signal

    fig, axs = plt.subplots(2, 2, figsize=(12, 10))

    # overview
    im = axs[0, 0].imshow(integrity_map, alpha=alpha, **plot_kwargs)
    boundaries_aligned.plot(ax=axs[0, 0], linewidth=0.1, **boundary_kwargs)
    axs[0, 0].set_title(sample_name)

    for (x, y, w, h) in boxes:
        rect = patches.Rectangle((x, y), w, h, linewidth=1, edgecolor='red',
                                 facecolor='none', linestyle='--')
        axs[0, 0].add_patch(rect)

    cax = make_axes_locatable(axs[0, 0]).append_axes("right", size="4%", pad=0.05)
    fig.colorbar(im, cax=cax).ax.set_title("VSI", fontsize=10)

    # histogram
    axs[0, 1].hist(vsi_mean, bins=100, zorder=3, color="slategrey")
    axs[0, 1].axvline(0.5, color='gray', linestyle='--', linewidth=1)
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
