import gzip
import shutil
import warnings
from typing import Dict

import dask.dataframe as dd
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio.features as rf
from dask.dataframe import DataFrame
from joblib import Parallel, delayed
from rasterio import Affine
from scipy import ndimage as ndi
from scipy.spatial import KDTree
from shapely.geometry import shape
from skimage.morphology import remove_small_holes, remove_small_objects
from skimage.segmentation import watershed
from tqdm import tqdm

# NOTE: spatialdata is imported lazily inside the functions that need it. Importing
# it at module top triggers spatialdata's dask-expr guard, which crashes joblib/loky
# workers (they re-import this module for _split_factor but never set the dask config).
from cellseg_benchmark import _constants

_F2CT = _constants.ficture_factor_to_celltype
_TRUE = _constants.true_cluster
_COLORS = _constants.cell_type_colors


def _resolve_cell_type(factor: pd.Series) -> pd.Series:
    """Map factor ids to their color-table cell-type names (see _constants.true_cluster)."""
    resolve = lambda ct: (_TRUE.get(ct) or [ct])[0]
    return factor.astype(str).map(_F2CT).map(resolve)


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


def load_pixel_tsv(
    file_path: str, skiprows: int = 3, chunksize: int = 10000
) -> pd.DataFrame:
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
    from spatialdata.models import PointsModel

    metadata = parse_metadata(pixel_level_factors_file)
    df = load_pixel_tsv(pixel_level_factors_file, skiprows=3)
    dask_df = dd.from_pandas(process_coordinates(df, metadata), npartitions=96)
    return PointsModel.parse(dask_df)


def get_transcript_level_factors(
    transcripts: pd.DataFrame,
    tree: KDTree,
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    current_factor: int,
) -> pd.DataFrame:
    """Assigns factor values to transcripts based on their nearest spatial location."""
    # query tree to get nearest pixels and according factor assignment
    query = np.array([transcripts["x"], transcripts["y"]]).T
    dd, ii = tree.query(query)
    # get factor prediction from df
    factor = np.array(df.iloc[ii]["K1"])
    # where distance > 5 um set factor to max_factor to indicate that this transcript was not mapped
    factor[dd > 5] = int(metadata["K"])
    return transcripts.assign(**{f"{current_factor}_factors": factor})


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


def build_factor_raster(pixel_file: str, res: float = 1.5, min_um2: float = 5.0):
    """Rasterize a FICTURE pixel decode to a majority-factor grid.

    Each grid cell stores its dominant factor + 1 (0 = background); tiny blobs
    are dropped and tiny holes filled.

    Args:
        pixel_file: Path to decode.pixel.sorted.tsv.gz.
        res: Grid size in micrometers.
        min_um2: Minimum blob / maximum hole area in um^2.

    Returns:
        tuple: (label raster, affine px->um, (xmin, xmax, ymin, ymax) in um).
    """
    meta = parse_metadata(pixel_file)
    offx, offy, scale = float(meta["OFFSET_X"]), float(meta["OFFSET_Y"]), float(meta["SCALE"])
    with gzip.open(pixel_file, "rt") as f:
        cols = next(l for l in f if l[:1] == "#" and l[:2] != "##")[1:].strip().split("\t")

    df = pd.read_csv(pixel_file, sep="\t", comment="#", names=cols, usecols=["X", "Y", "K1"])
    x, y = df.X / scale + offx, df.Y / scale + offy
    grid = pd.DataFrame({"row": ((y - offy) / res).astype(int),
                         "col": ((x - offx) / res).astype(int), "K1": df.K1})
    counts = grid.groupby(["row", "col", "K1"]).size().reset_index(name="n")
    win = counts.loc[counts.groupby(["row", "col"])["n"].idxmax()]
    lab = np.zeros((int(win.row.max()) + 1, int(win.col.max()) + 1), np.int32)
    lab[win.row, win.col] = win.K1 + 1

    mp = max(1, int(min_um2 / res ** 2))
    clean = np.zeros_like(lab)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        for fac in np.unique(lab[lab > 0]):
            m = remove_small_holes(remove_small_objects(lab == fac, mp), area_threshold=mp)
            clean[m & (clean == 0)] = fac

    affine = Affine(res, 0, offx, 0, res, offy)
    bbox = (float(x.min()), float(x.max()), float(y.min()), float(y.max()))
    return clean, affine, bbox


def segments_to_boundaries(lab: np.ndarray, affine: Affine) -> gpd.GeoDataFrame:
    """Polygonize raw FICTURE segments (one polygon per connected same-factor patch).

    Args:
        lab: Majority-factor raster from build_factor_raster.
        affine: Grid->um transform.

    Returns:
        GeoDataFrame indexed by segment_id with factor, cell_type, area_um2.
    """
    gdf = gpd.GeoDataFrame(
        [{"factor": int(v) - 1, "geometry": shape(g)}
         for g, v in rf.shapes(lab, mask=lab > 0, transform=affine, connectivity=8)])
    gdf["cell_type"] = _resolve_cell_type(gdf.factor)
    gdf["area_um2"] = gdf.geometry.area
    gdf.index.name = "segment_id"
    return gdf


def _split_factor(sub, rr, cc, ee, conn):
    """Split one factor crop into cells (module-level so joblib/loky can pickle it).

    Each mask pixel is assigned to its nearest nucleus within the crop (geodesic
    Voronoi = seeded watershed on a flat image); nucleus-free components are kept.

    Returns:
        tuple: (label crop with local ids, [(entity_id, n_nuclei) per local id]).
    """
    out = np.zeros(sub.shape, np.int32)
    metas, lab_ws = [], np.zeros(sub.shape, np.int32)
    if len(rr):
        markers = np.zeros(sub.shape, np.int32)
        for i, (r, c) in enumerate(zip(rr, cc), start=1):
            markers[r, c] = i
            metas.append((ee[i - 1], 1))
        lab_ws = watershed(np.zeros(sub.shape, np.uint8), markers=markers,
                           mask=sub, connectivity=conn)
        out[lab_ws > 0] = lab_ws[lab_ws > 0]
    nid = len(metas) + 1
    comp, ncomp = ndi.label(sub & (lab_ws == 0))
    for k in range(1, ncomp + 1):
        out[comp == k] = nid
        metas.append((None, 0))
        nid += 1
    return out, metas


def split_by_nuclei(lab: np.ndarray, affine: Affine, nuclei_xy, entity_ids,
                    connectivity: int = 1, n_jobs: int = 8) -> gpd.GeoDataFrame:
    """Split FICTURE segments into single cells using DAPI nuclei as seeds.

    Each pixel is assigned to its nearest nucleus within its factor segment
    (geodesic Voronoi); segments without a nucleus are kept and flagged
    n_nuclei == 0. Factors are processed in parallel.

    Args:
        lab: Majority-factor raster from build_factor_raster.
        affine: Grid->um transform (a/c/f give res, offset_x, offset_y).
        nuclei_xy: (N, 2) nucleus centroids in um.
        entity_ids: (N,) nucleus ids aligned to nuclei_xy.
        connectivity: 1 (4-conn, avoids diagonal 1px bridges) or 2 (8-conn).
        n_jobs: Parallel workers.

    Returns:
        GeoDataFrame indexed by cell_id with factor, entity_id, n_nuclei,
        cell_type, area_um2.
    """
    res, offx, offy = affine.a, affine.c, affine.f
    xy = np.asarray(nuclei_xy, float)
    rows = ((xy[:, 1] - offy) / res).astype(int)
    cols = ((xy[:, 0] - offx) / res).astype(int)
    inb = (rows >= 0) & (rows < lab.shape[0]) & (cols >= 0) & (cols < lab.shape[1])
    if inb.mean() < 0.5:  # wrong frame -> most nuclei fall outside the raster
        raise ValueError(
            f"nuclei and pixel coordinate frames likely differ: only "
            f"{inb.mean():.0%} of nuclei fall within the raster")
    rows, cols, ids = rows[inb], cols[inb], np.asarray(entity_ids)[inb]

    boxes = ndi.find_objects(lab)
    nuc_fac = lab[rows, cols]
    tasks = [(int(fac), sl := boxes[fac - 1], lab[sl] == fac,
              rows[nuc_fac == fac] - sl[0].start, cols[nuc_fac == fac] - sl[1].start,
              ids[nuc_fac == fac]) for fac in np.unique(lab[lab > 0])]
    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_split_factor)(sub, rr, cc, ee, connectivity)
        for _, _, sub, rr, cc, ee in tasks)

    ws = np.zeros_like(lab)
    factor, entity, nnuc = {}, {}, {}
    cid = 0
    for (fac, sl, *_), (out, metas) in zip(tasks, results):
        ws[sl][out > 0] = out[out > 0] + cid   # offset local ids -> global
        for j, (e, nn) in enumerate(metas, start=1):
            factor[cid + j], entity[cid + j], nnuc[cid + j] = fac - 1, e, nn
        cid += len(metas)

    gdf = gpd.GeoDataFrame(
        [{"cell_id": int(v), "factor": factor[int(v)], "entity_id": entity[int(v)],
          "n_nuclei": nnuc[int(v)], "geometry": shape(g)}
         for g, v in rf.shapes(ws, mask=ws > 0, transform=affine, connectivity=8)]
    ).dissolve(by="cell_id", aggfunc={"factor": "first", "entity_id": "first",
                                      "n_nuclei": "first"})
    gdf["cell_type"] = _resolve_cell_type(gdf.factor)
    gdf["area_um2"] = gdf.geometry.area
    return gdf


def write_boundaries(gdf: gpd.GeoDataFrame, zarr) -> None:
    """Write polygons as shapes["boundaries"] in a fresh SpatialData .zarr (Identity/micron)."""
    import spatialdata as sd
    from spatialdata.models import ShapesModel
    from spatialdata.transformations import Identity

    shapes = ShapesModel.parse(gdf, transformations={"micron": Identity()})
    sd.SpatialData(shapes={"boundaries": shapes}).write(str(zarr), overwrite=True)


def aggregate_tables(data_path: str, targets, gene_column: str = "gene",
                     min_transcripts: int = 10, n_workers: int = 4) -> None:
    """Assign MERSCOPE transcripts to boundary sets and write each as sdata with a table.

    Loads the transcripts once (sopa.io.merscope) and aggregates counts per
    polygon for every target (as in the repo's segmentation scripts). Each
    output holds shapes["boundaries"] and the cell x gene AnnData under
    tables["table"]; polygons below min_transcripts are dropped (shapes and
    table stay consistent).

    Args:
        data_path: MERSCOPE output folder.
        targets: iterable of (GeoDataFrame in um, output sdata.zarr path).
        gene_column: Transcript gene column name.
        min_transcripts: Minimum transcripts per polygon to keep.
        n_workers: Dask workers for aggregation.
    """
    import anndata._core.anndata as anndata_core
    import scipy.sparse as sp
    import sopa
    import spatialdata as sd
    from spatialdata.models import ShapesModel
    from spatialdata.transformations import Identity, get_transformation

    # sopa 2.0.6 builds COO count matrices; anndata >=0.12 only accepts CSR/CSC
    _coerce = anndata_core.coerce_array
    anndata_core.coerce_array = lambda v, **kw: _coerce(
        v.tocsr() if sp.issparse(v) and v.format not in ("csr", "csc") else v, **kw)

    src = sopa.io.merscope(data_path)
    tx = src[list(src.points.keys())[0]]
    cs = list(get_transformation(tx, get_all=True).keys())

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = n_workers
    for gdf, zarr in targets:
        sdata = sd.SpatialData(
            points={"transcripts": tx},
            shapes={"boundaries": ShapesModel.parse(
                gdf, transformations={c: Identity() for c in cs})},
        )
        sdata.attrs["transcripts_dataframe"] = "transcripts"

        tmp = f"{zarr}.tmp"                                 # sopa needs a backed store
        sdata.write(tmp, overwrite=True)
        sdata = sd.read_zarr(tmp)
        sopa.aggregate(sdata, gene_column=gene_column, aggregate_channels=False,
                       min_transcripts=min_transcripts, shapes_key="boundaries")
        del sdata["transcripts"]
        sdata.write(str(zarr), overwrite=True)
        shutil.rmtree(tmp, ignore_errors=True)


def plot_qc(gdf: gpd.GeoDataFrame, nuclei, aspect: float, title: str, path) -> None:
    """Save a full-slide QC image: cells in cell-type colors, nuclei as short crosses."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(30, 30 * aspect))
    gdf.plot(ax=ax, color=gdf.cell_type.map(_COLORS), edgecolor="black", linewidth=0.1)
    if nuclei is not None:
        ax.plot(nuclei[:, 0], nuclei[:, 1], "+", color="black",
                markersize=0.9, markeredgewidth=0.2)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title)
    fig.savefig(str(path), dpi=300, bbox_inches="tight")
    plt.close(fig)
