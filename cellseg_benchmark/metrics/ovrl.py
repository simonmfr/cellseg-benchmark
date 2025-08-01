import logging
import os
from os.path import exists, join
from typing import List, Tuple

import pandas as pd
from dask.array import from_array
from geopandas import GeoDataFrame
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
