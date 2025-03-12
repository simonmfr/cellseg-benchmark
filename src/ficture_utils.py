import gzip
import pandas as pd
import dask.dataframe as dd
from tqdm import tqdm
from spatialdata.models import PointsModel


def parse_metadata(file_path):
    """
    Parse metadata from the first three lines of FICTURE pixel-level tsv.gz file.

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


def load_pixel_tsv(file_path, skiprows=3, chunksize=10000):
    """
    Load FICTURE pixel-level tsv.gz file with progress bar.

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


def process_coordinates(df, metadata):
    """
    Transform FICTURE pixel coordinates to micrometer scale and rename columns.

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


def get_pixel_level_factors(pixel_level_factors_file):
    """
    Load and format FICTURE pixel-level file to micrometer scale and return a parsed SpatialData PointsModel.

    Args:
        pixel_level_factors_file (str): Path to the FICTURE pixel-level tsv.gz file (*.pixel.sorted.tsv.gz)

    Returns:
        PointsModel: Dask dataframe parsed  as a SpatialData PointsModel.
    """
    metadata = parse_metadata(pixel_level_factors_file)
    df = load_pixel_tsv(pixel_level_factors_file, skiprows=3)
    dask_df = dd.from_pandas(process_coordinates(df, metadata), npartitions=96)
    return PointsModel.parse(dask_df)


def get_transcript_level_factors(transcripts, tree, df, metadata, current_factor):
    # query tree to get nearest pixels and according factor assignment
    query = np.array([transcripts["x"], transcripts["y"]]).T
    dd, ii = tree.query(query)
    # get factor prediction from df
    factor = np.array(df.iloc[ii]["K1"])
    # where distance > 5 um set factor to max_factor to indicate that this transcript was not mapped
    factor[dd > 5] = int(metadata["K"])
    return transcripts.assign(**{f"{current_factor}_factors": factor})
