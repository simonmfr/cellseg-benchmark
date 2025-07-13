import gzip
from typing import Dict

from dask.dataframe import DataFrame
import dask.dataframe as dd
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from spatialdata.models import PointsModel
from tqdm import tqdm


def parse_metadata(file_path: str) -> Dict[str, str]:
    """Parse metadata from the first three lines of FICTURE pixel-level tsv.gz file.

    Args:
        file_path (str): Path to file.

    Returns:
        dict: Extracted metadata as key-value pairs.
    """
    metadata = {}
    with gzip.open(file_path, "rb") as f:
        for i, line in enumerate(f):
            if i >= 3:
                break
            line = line.decode().strip("#").strip()
            for s in line.split(";"):
                k, v = s.split("=")
                metadata[k] = v
    return metadata


def load_pixel_tsv(file_path: str, skiprows: int=3, chunksize: int=10000) -> pd.DataFrame:
    """Load FICTURE pixel-level tsv.gz file with progress bar.

    Args:
        file_path (str): Path to the FICTURE pixel-level tsv.gz file (*.pixel.sorted.tsv.gz)
        skiprows (int): Number of rows to skip at the beginning of the file.
        chunksize (int): Number of rows to read per chunk.

    Returns:
        pd.DataFrame: Concatenated DataFrame from all chunks.
    """
    with gzip.open(file_path, "rb") as f:
        total_lines = sum(1 for _ in f)
        f.seek(0)
        df_chunks = pd.read_table(
            f,
            skiprows=skiprows,
            header=0,
            engine="c",
            iterator=True,
            chunksize=chunksize,
        )
        return pd.concat(
            tqdm(df_chunks, total=total_lines // chunksize, desc="Loading data"),
            ignore_index=True,
        )


def process_coordinates(df: pd.DataFrame, metadata: Dict[str, str]) -> pd.DataFrame:
    """Transform FICTURE pixel coordinates to micrometer scale and rename columns.

    Args:
        df (pd.DataFrame): DataFrame containing pixel coordinates.
        metadata (dict): Metadata containing scale and offset values.

    Returns:
        pd.DataFrame: Updated DataFrame with transformed and renamed coordinates.
    """
    scale = float(metadata["SCALE"])
    offset_x = float(metadata["OFFSET_X"])
    offset_y = float(metadata["OFFSET_Y"])

    df["X_um"] = df["X"] / scale + offset_x
    df["Y_um"] = df["Y"] / scale + offset_y
    df = df.sort_values(["X_um", "Y_um"])
    return df.rename(columns={"X_um": "x", "Y_um": "y", "X": "X_px", "Y": "Y_px"})


def get_pixel_level_factors(pixel_level_factors_file: str) -> DataFrame:
    """Load and format FICTURE pixel-level file to micrometer scale and return a parsed SpatialData PointsModel.

    Args:
        pixel_level_factors_file (str): Path to the FICTURE pixel-level tsv.gz file (*.pixel.sorted.tsv.gz)

    Returns:
        PointsModel: Dask dataframe parsed  as a SpatialData PointsModel.
    """
    metadata = parse_metadata(pixel_level_factors_file)
    df = load_pixel_tsv(pixel_level_factors_file, skiprows=3)
    dask_df = dd.from_pandas(process_coordinates(df, metadata), npartitions=96)
    return PointsModel.parse(dask_df)


def get_transcript_level_factors(transcripts: pd.DataFrame, tree: KDTree, df: pd.DataFrame, metadata: pd.DataFrame, current_factor: int) -> pd.DataFrame:
    """Assigns factor values to transcripts based on their nearest spatial location."""
    # query tree to get nearest pixels and according factor assignment
    query = np.array([transcripts["x"], transcripts["y"]]).T
    dd, ii = tree.query(query)
    # get factor prediction from df
    factor = np.array(df.iloc[ii]["K1"])
    # where distance > 5 um set factor to max_factor to indicate that this transcript was not mapped
    factor[dd > 5] = int(metadata["K"])
    return transcripts.assign(**{f"{current_factor}_factors": factor})


def create_factor_level_image(data, factor, DAPI_shape, top_n_factors: int) -> np.ndarray:
    """Compute image for given factor.

    Args:
        data: ficture data
        factor: factor to compute image for
        DAPI_shape: target shape

    Returns: image of factor

    """
    filtered_data = []
    for i in range(1, top_n_factors + 1):
        filtered_data.append(pd.concat([
            data.loc[data[f"K{i}"] == factor, ["Y_pixel", "X_pixel"]],
            data.loc[data[f"K{i}"] == factor, f"P{i}"].rename(f"P{i}": "probability"
                                                              , inplace=False)
        ], axis=1)
        )
    filtered_data = pd.concat(filtered_data, axis=0)

#    K2_ind = data["K2"] == factor
#    K2 = data[K2_ind]
#    K2["probability"] = K2["P2"].copy()

#    K3_ind = data["K3"] == factor
#    K3 = data[K3_ind]
#    K3["probability"] = K3["P3"].copy()

#    filtered_data = pd.concat([K1, K2, K3], axis=0)[
#        ["Y_pixel", "X_pixel", "probability"]
#    ]
#    del K1, K2, K3
#    filtered_data = K1[["Y_pixel", "X_pixel", "probability"]]
#    del K1

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
