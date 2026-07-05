#!/usr/bin/env python
"""FICTURE output I/O shared across the package (metrics and master-sdata builds)."""

import gzip
import pathlib
from typing import Dict

import dask
import numpy as np
import pandas as pd

dask.config.set(
    {"dataframe.query-planning": False}
)  # prevents a seg_postprocessing runtime error


def parse_metadata(file_path: str) -> Dict[str, str]:
    """Parse the ``key=value`` metadata from the first three lines of a FICTURE pixel file."""
    metadata = {}
    with gzip.open(file_path, "rb") as f:
        for i, line in enumerate(f):
            if i >= 3:
                break
            for entry in line.decode().strip("#").strip().split(";"):
                key, value = entry.split("=")
                metadata[key] = value
    return metadata


def read_ficture_pixels(
    path, header=("BLOCK", "X", "Y", "K1", "K2", "K3", "P1", "P2", "P3"), usecols=None
) -> pd.DataFrame:
    """Read a FICTURE pixel-level decode file (``decode.pixel.sorted.tsv.gz``).

    Pass ``usecols`` (e.g. ``["X", "Y", "K1"]``) to load only the needed columns.
    """
    return pd.read_csv(
        path, sep="\t", names=list(header), comment="#", usecols=usecols
    )


def find_ficture_output(sample: str, base_path: str) -> str:
    """Return the path to a sample's FICTURE decode pixel file.

    Since the FICTURE pipeline refactor the output layout is fixed at
    ``samples/<sample>/results/Ficture/output/decode.pixel.sorted.tsv.gz``.

    Args:
        sample: Sample name.
        base_path: Benchmark base directory.

    Returns:
        Path to ``decode.pixel.sorted.tsv.gz``.

    Raises:
        FileNotFoundError: If the file does not exist for the sample.
    """
    path = (
        pathlib.Path(base_path)
        / "samples"
        / sample
        / "results"
        / "Ficture"
        / "output"
        / "decode.pixel.sorted.tsv.gz"
    )
    if not path.exists():
        raise FileNotFoundError(f"No FICTURE output for sample '{sample}': {path}")
    return str(path)


def create_factor_level_image(
    data, factor, DAPI_shape, top_n_factors: int
) -> np.ndarray:
    """Compute image for given factor.

    Args:
        data: ficture data
        factor: factor to compute image for
        DAPI_shape: target shape
        top_n_factors: number of top factors work with

    Returns: image of factor

    """
    filtered_data = []
    for i in range(1, top_n_factors + 1):
        filtered_data.append(
            pd.concat(
                [
                    data.loc[data[f"K{i}"] == factor, ["Y_pixel", "X_pixel"]],
                    data.loc[data[f"K{i}"] == factor, f"P{i}"].rename(
                        "probability", inplace=False
                    ),
                ],
                axis=1,
            )
        )
    filtered_data = pd.concat(filtered_data, axis=0)

    bins_y = np.linspace(0, DAPI_shape[1], num=DAPI_shape[1] + 1)
    bins_x = np.linspace(0, DAPI_shape[0], num=DAPI_shape[0] + 1)
    image, _, _ = np.histogram2d(
        filtered_data["Y_pixel"],
        filtered_data["X_pixel"],
        bins=[bins_x, bins_y],
        weights=filtered_data["probability"],
    )
    image = np.clip(
        np.around(
            image * (np.finfo(np.float16).max.astype(np.uint16) - 5)
        ),  # ensures no overflow with np.float16
        0,
        (np.finfo(np.float16).max.astype(np.uint16) - 5),
    ).astype(np.uint16)  # makes smaller file
    image = image[np.newaxis, :]
    return image
