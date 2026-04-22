from pathlib import Path

import dask
import geopandas as gpd
import numpy as np
import numpy.ma as ma
import pandas as pd
import shapely
import spatialdata as sd
from anndata import AnnData
from dask.diagnostics import ProgressBar
from joblib import delayed, Parallel
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, box
from sopa.segmentation.shapes import pixel_outer_bounds, rasterize
from xarray import DataArray

from cellseg_benchmark._constants import true_cluster, factor_to_celltype
from cellseg_benchmark.ficture_utils import (
    process_coordinates,
    parse_metadata,
    read_ficture_pixels,
    find_ficture_output
)
from cellseg_benchmark.metrics.f1_score import compute_f1
from cellseg_benchmark.metrics.utils import _prepare_label_maps
from cellseg_benchmark.sdata_utils import _assign_points_to_polygons

def _process_sample_ficture_f1(
    sample,
    obs_df,
    method,
    base_path,
    n_ficture,
    factor_to_celltype,
    true_cluster,
):
    """Compute ficture stats on one sample."""
    ficture_full_path = find_ficture_output(sample, base_path, n_ficture)

    sdata = sd.read_zarr(
        Path(base_path) / "samples" / sample / "sdata_z3.zarr",
        selection=("points", "shapes"),
    )

    points_data = sdata[f"{sample}_transcripts"].compute()
    boundaries = sd.transform(
        sdata[f"boundaries_{method}"],
        to_coordinate_system="micron",
    )

    del sdata
    boundaries = boundaries.copy()
    boundaries["p_id"] = boundaries.index

    ficture_pixels = read_ficture_pixels(ficture_full_path)

    metadata = parse_metadata(ficture_full_path)
    ficture_pixels = process_coordinates(ficture_pixels, metadata)

    fic_microns = ficture_pixels[["x", "y"]].to_numpy()
    points_coords = points_data[["x", "y"]].to_numpy()

    tree = cKDTree(fic_microns)
    distances, indices = tree.query(points_coords, k=1, distance_upper_bound=5)

    result = points_data.copy()
    valid = np.isfinite(distances)
    result["assigned_factor"] = -1
    result["nearest_distance"] = distances
    result.loc[valid, "assigned_factor"] = ficture_pixels.iloc[indices[valid]][
        "K1"
    ].to_numpy()
    del ficture_pixels

    res = _assign_points_to_polygons(
        result,
        boundaries,
        polygon_id_col="p_id",
        output_col="assigned_polygon",
    )

    celltype = obs_df.loc[obs_df["sample"] == sample, "cell_type_revised"].copy()

    res_cts = res.merge(
        celltype,
        how="left",
        left_on="assigned_polygon",
        right_index=True,
    )

    factor_map, correct_celltypes = _prepare_label_maps(
        factor_to_celltype, true_cluster
    )

    res_cts["assigned_factor"] = res_cts["assigned_factor"].astype(str).map(factor_map)
    res_cts["cell_type_revised"] = res_cts["cell_type_revised"].map(correct_celltypes)
    eval_df = res_cts[["cell_type_revised", "assigned_factor"]].copy()

    metric = compute_f1(
        eval_df,
        celltype_col="cell_type_revised",
        factor_col="assigned_factor",
        flavor="all",
        return_confusion=True,
    )

    sample_results = metric.T.reset_index().rename(
        columns={
            "index": "cell_type",
            "F1_score": "f1_score",
            "TP": "tp",
            "FP": "fp",
            "FN": "fn",
        }
    )
    sample_results = sample_results[
        ~sample_results["cell_type"].isin(["macro F1_score", "micro F1_score"])
    ].copy()

    sample_results["precision"] = sample_results["tp"] / (
        sample_results["tp"] + sample_results["fp"]
    )
    sample_results["recall"] = sample_results["tp"] / (
        sample_results["tp"] + sample_results["fn"]
    )

    sample_results["precision"] = sample_results["precision"].fillna(0)
    sample_results["recall"] = sample_results["recall"].fillna(0)
    sample_results.insert(0, "sample", sample)

    return sample, sample_results, eval_df


def ficture_f1_parallel(
        adata: AnnData,
        method: str,
        base_path: str | Path,
        n_ficture: int,
        n_jobs: int = -1,
        backend: str = "loky",
        sample_col: str = "sample",
        celltype_col: str = "cell_type_revised",
) -> pd.DataFrame:
    """
    Compute Ficture F1 score with parallelization.

    Args:
        adata (AnnData): adata containing sample and cell type columns in .obs.
        method (str): method to compute ficture F1 score for on all samples in adata.obs.
        base_path (str): Base data directory.
        n_ficture (int): number of fictures factors.
        n_jobs (int): number of parallel jobs.
        backend (str): backend of parallelization.
        sample_col (str): column in adata.obs containing sample names.
        celltype_col (str): column in adata.obs containing cell type names.

    Returns:
        DataFrame with columns [...] #TODO
    """

    obs_df = adata.obs[[sample_col, celltype_col]].copy()
    obs_df.index = [x[:10] for x in adata.obs_names]

    samples = obs_df[sample_col].unique().tolist()

    worker_out = Parallel(n_jobs=n_jobs, backend=backend, verbose=10)(
        delayed(_process_sample_ficture_f1)(
            sample=sample,
            obs_df=obs_df,
            method=method,
            base_path=base_path,
            n_ficture=n_ficture,
            factor_to_celltype=factor_to_celltype,
            true_cluster=true_cluster,
        )
        for sample in samples
    )

    per_sample_tables = []
    eval_frames = {}

    for sample, sample_results, eval_df in worker_out:
        per_sample_tables.append(sample_results)
        eval_frames[sample] = eval_df

    results = pd.concat(per_sample_tables, ignore_index=True)

    # Global score across all samples
    overall_metric = compute_f1(
        eval_frames,
        celltype_col=celltype_col,
        factor_col="assigned_factor",
        flavor="all",
        return_confusion=True,
    )

    overall_results = overall_metric.T.reset_index().rename(
        columns={
            "index": "cell_type",
            "F1_score": "f1_score",
            "TP": "tp",
            "FP": "fp",
            "FN": "fn",
        }
    )
    overall_results = overall_results[
        ~overall_results["cell_type"].isin(["macro F1_score", "micro F1_score"])
    ].copy()

    overall_results["precision"] = overall_results["tp"] / (
            overall_results["tp"] + overall_results["fp"]
    )
    overall_results["recall"] = overall_results["tp"] / (
            overall_results["tp"] + overall_results["fn"]
    )

    overall_results["precision"] = overall_results["precision"].fillna(0)
    overall_results["recall"] = overall_results["recall"].fillna(0)
    overall_results.insert(0, "sample", "all_samples")

    results = pd.concat([results, overall_results], ignore_index=True)
    results.drop(results[results["cell_type"] == "unassigned"].index, inplace=True)

    return results

# from sopa.aggregation.channels.py
AVAILABLE_MODES = ["average", "min", "max", "variance", "sum"]

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
