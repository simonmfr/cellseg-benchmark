from __future__ import annotations

from pathlib import Path
from typing import Iterable, Dict

import dask
import pandas as pd
from dask.diagnostics import ProgressBar
import geopandas as gpd
import numpy as np
import numpy.ma as ma
from shapely.affinity import translate, scale
from shapely.geometry import box
from shapely import affinity, STRtree
from dask.base import is_dask_collection
import numpy as np
import dask


# --- SpatialData stack (adapt if your imports differ)
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Affine, set_transformation
from spatialdata import read_zarr
from sopa.utils import get_spatial_image, to_intrinsic  # <- you already have these
from sopa.segmentation.shapes import expand_radius, pixel_outer_bounds, rasterize

from ._constants import image_based
from .sdata_utils import prepare_ficture

def _to_shapely_coeffs(M: np.ndarray) -> list[float]:
    """Convert 3x3 homogeneous affine to shapely [a,b,d,e,xoff,yoff]."""
    return [float(M[0,0]), float(M[0,1]), float(M[1,0]), float(M[1,1]), float(M[0,2]), float(M[1,2])]

def build_boundary_gdf_in_micron(
    method: str,
    zarr_path: Path,
    transformation: np.ndarray,
    boundary_key: str | None,
) -> tuple[str, gpd.GeoDataFrame]:
    """Load shapes for one method and return a GDF in micron coords (no sdata mutation; safe for parallel)."""
    tmp = read_zarr(zarr_path, selection=("shapes", "tables"))
    if boundary_key is None:
        boundary_key = tmp[list(tmp.tables.keys())[0]].uns["spatialdata_attrs"]["region"]

    # mirror your vpt_3D dissolving behavior
    if method.startswith("vpt_3D"):
        try:
            gdf = tmp[boundary_key][["cell_id", "geometry"]].dissolve(by="cell_id")
        except ValueError:
            df = tmp[boundary_key][["cell_id", "geometry"]].copy()
            df.index.rename(None, inplace=True)
            gdf = df.dissolve(by="cell_id")
    else:
        if "cell_id" in tmp[boundary_key].columns:
            gdf = tmp[boundary_key][["cell_id", "geometry"]].copy()
        else:
            gdf = tmp[boundary_key][["geometry"]].copy()
            gdf["cell_id"] = gdf.index

    # put geometries into MICRON (matches your set_transformation logic)
    needs_inverse = any(method.startswith(p) for p in image_based) and method != "Cellpose_1_Merlin"
    if needs_inverse:
        Minv = np.linalg.inv(transformation)
        coeffs = _to_shapely_coeffs(Minv)
        gdf["geometry"] = gdf.geometry.apply(lambda geom: affinity.affine_transform(geom, coeffs))
    # else: Identity to micron (no-op)

    if "cell_id" in gdf.columns:
        gdf = gdf.set_index("cell_id", drop=True)

    # ensure GeoDataFrame
    gdf = gpd.GeoDataFrame(gdf, geometry="geometry", crs=None)
    return method, gdf

def build_boundaries_parallel(
    methods: Iterable[str],
    results_path: Path,
    transformation: np.ndarray,
    boundary_key: str | None = None,
    scheduler: str | None = None,
) -> Dict[str, gpd.GeoDataFrame]:
    """Return dict(method -> GDF in micron) by parallelizing loading/dissolving/transforming across methods."""
    tasks = [
        dask.delayed(build_boundary_gdf_in_micron)(
            m, results_path / m / "sdata.zarr", transformation, boundary_key
        )
        for m in methods
    ]
    with ProgressBar():
        out = dask.compute(*tasks, scheduler=scheduler)
    return dict(out)

def _tiles_and_intersections(cells, chunk_sizes_y, chunk_sizes_x):
    """Compute (on driver) which cells intersect which image tiles."""
    offsets_y = np.cumsum(np.pad(chunk_sizes_y, (1, 0)))
    offsets_x = np.cumsum(np.pad(chunk_sizes_x, (1, 0)))

    tree = STRtree(cells)
    tiles = []
    for iy in range(len(chunk_sizes_y)):
        ymin, ymax = offsets_y[iy], offsets_y[iy + 1]
        for ix in range(len(chunk_sizes_x)):
            xmin, xmax = offsets_x[ix], offsets_x[ix + 1]
            patch = box(xmin, ymin, xmax, ymax)
            idxs = tree.query(patch, predicate="intersects")
            if len(idxs): #only add patch when overlaps exist
                tiles.append((iy, ix, (xmin, xmax, ymin, ymax), idxs))
    return tiles

def _aggregate_chunk_task(
    chunk: np.ndarray,           # (C, h, w)
    tile_bounds,                 # (xmin, xmax, ymin, ymax)
    tile_idxs: list[int],        # indices into global cells list
    cells_subset,                # polygons for those indices
    mode: str,
    means_subset: np.ndarray | None,  # (n_local, C) for variance
):
    C = chunk.shape[0]
    n_local = len(tile_idxs)
    if mode == "min":
        local_agg = np.full((n_local, C), np.inf, dtype=np.float64)
    elif mode == "max":
        local_agg = np.full((n_local, C), -np.inf, dtype=np.float64)
    else:
        local_agg = np.zeros((n_local, C), dtype=np.float64)
    local_area = np.zeros(n_local, dtype=np.float64)

    xmin, xmax, ymin, ymax = tile_bounds

    for j, cell in enumerate(cells_subset):
        cb = pixel_outer_bounds(cell.bounds)  # (xmin,ymin,xmax,ymax) in global pixel coords

        y0 = max(cb[1] - ymin, 0); y1 = min(cb[3] - ymin, chunk.shape[1])
        x0 = max(cb[0] - xmin, 0); x1 = min(cb[2] - xmin, chunk.shape[2])
        if y1 <= y0 or x1 <= x0:
            continue

        sub = chunk[:, y0:y1, x0:x1]
        mask = rasterize(cell, sub.shape[1:], cb)
        if mask.size == 0:
            continue

        local_area[j] += mask.sum()

        if mode == "min":
            masked = ma.masked_array(sub, 1 - np.repeat(mask[None], C, axis=0))
            local_agg[j] = np.minimum(local_agg[j], masked.min(axis=(1, 2)).filled(np.inf).astype(np.float64))
        elif mode == "max":
            local_agg[j] = np.maximum(local_agg[j], (sub * mask).max(axis=(1, 2)).astype(np.float64))
        elif mode in ("sum", "average"):
            vals = np.sum(sub * mask, axis=(1, 2), dtype=np.uint64)
            local_agg[j] += vals
        elif mode == "variance":
            # means_subset[j] shape: (C,)
            sub_f = sub.astype(np.float32) * (1.0 / 65535.0)
            diff = sub_f - means_subset[j][:, None, None]  # (C, h, w)
            vals = (diff * diff * mask).sum(axis=(1, 2), dtype=np.float64)
            local_agg[j] += vals

    return tile_idxs, local_agg, local_area


def _aggregate_channels_aligned_parallel_lazy(
    image_da,                        # xarray.DataArray with dask backend, shape (c, y, x)
    geo_df,                          # GeoDataFrame; geometries already in intrinsic (pixel) space
    mode: str,
    means_delayed=None,              # np.ndarray | Future | delayed (n_cells, C) if variance
):
    """
    Parallel, tile-wise aggregation with late scaling.
    - If a distributed Client is active, returns a distributed.Future.
    - Otherwise returns a dask.delayed object.

    Avoids 'bytes object is too large' by scattering big, passing small:
      * scatter the full list of cells (polygons) once per worker
      * pass only 'idxs' (ints) to each tile task and slice means per tile on the worker
    """
    assert mode in {"average", "max", "min", "variance", "sum"}

    # --- derive dims/axes robustly, and ensure one chunk along 'c' ---
    image_da = image_da.chunk({"c": -1})
    darr = image_da.data
    dims = list(image_da.dims)              # expect ["c","y","x"] but don't assume order
    idx = {d: i for i, d in enumerate(dims)}
    c_ax, y_ax, x_ax = idx["c"], idx["y"], idx["x"]

    # chunk sizes for tiling
    chunk_sizes_y = darr.chunks[y_ax]
    chunk_sizes_x = darr.chunks[x_ax]

    # tiles & blocks
    cells = list(geo_df.geometry)
    tiles = _tiles_and_intersections(cells, chunk_sizes_y, chunk_sizes_x)
    blocks = np.array(darr.to_delayed())    # shape = numblocks per axis

    # scale for late normalization
    dt = np.dtype(image_da.dtype)
    if np.issubdtype(dt, np.integer):
        scale = float(np.iinfo(dt).max) if np.iinfo(dt).max > 1 else 1.0
    else:
        scale = 1.0  # assume already normalized floats

    # --- reducer: combine per-tile partials into final array ---
    def _reduce_parts(parts, n_cells, n_channels, mode, scale):
        if mode == "min":
            agg = np.full((n_cells, n_channels), np.inf, dtype=np.float64)
        elif mode == "max":
            agg = np.full((n_cells, n_channels), -np.inf, dtype=np.float64)
        else:
            agg = np.zeros((n_cells, n_channels), dtype=np.float64)
        areas = np.zeros(n_cells, dtype=np.float64)

        for idxs, local_agg, local_area in parts:
            idxs = np.asarray(idxs)
            if mode == "max":
                np.maximum.at(agg, (idxs, slice(None)), local_agg)
            elif mode == "min":
                np.minimum.at(agg, (idxs, slice(None)), local_agg)
            else:
                np.add.at(agg, (idxs, slice(None)), local_agg)
            np.add.at(areas, idxs, local_area)

        if mode == "average":
            denom = np.clip(areas, 1, None)[:, None]
            return (agg / denom / scale).astype(np.float32)
        elif mode == "variance":
            # already accumulated sum((x_norm - mean_norm)^2) in tiles
            denom = np.clip(areas[:, None] - 1, 1, None)
            return (agg / denom).astype(np.float32)
        elif mode in ("max", "min"):
            return (agg / scale).astype(np.float32)
        elif mode == "sum":
            return agg.astype(np.float64)

    # --- worker-side wrapper: reconstruct small per-tile inputs from scattered big ones ---
    def _run_one_scattered(chunk, tile_bounds, idxs, mode, means_or_none, cells_all):
        # 'cells_all' and 'means_or_none' are realized on the worker (scattered or computed)
        cells_subset = [cells_all[k] for k in idxs]
        means_subset = None
        if mode == "variance" and means_or_none is not None:
            means_subset = means_or_none[np.asarray(idxs)]
        return _aggregate_chunk_task(
            chunk, tile_bounds, idxs, cells_subset, mode, means_subset
        )

    # try to use a distributed client; otherwise fall back to delayed path
    try:
        from dask.distributed import get_client
        client = get_client()
    except Exception:
        client = None

    if client is not None:
        # ---- distributed branch: scatter big, pass small, submit futures ----
        cells_fut = client.scatter(cells, broadcast=True)
        mean_ref = None
        if mode == "variance":
            if means_delayed is None:
                raise ValueError("means required for variance computation")
            if is_dask_collection(means_delayed):
                mean_ref = client.compute(means_delayed)     # compute once per method
            else:
                # numpy array or already a Future
                try:
                    from distributed.future import Future  # type: ignore
                    is_future = isinstance(means_delayed, Future)
                except Exception:
                    is_future = False
                mean_ref = means_delayed if is_future else client.scatter(means_delayed, broadcast=True)

        futures = []
        for (iy, ix, (xmin, xmax, ymin, ymax), idxs_) in tiles:
            bi = [0] * darr.ndim
            bi[c_ax] = 0            # channels fused to one block
            bi[y_ax] = iy
            bi[x_ax] = ix
            chunk_ref = blocks[tuple(bi)]
            fut = client.submit(
                _run_one_scattered,
                chunk_ref,
                (xmin, xmax, ymin, ymax),
                idxs_,
                mode,
                mean_ref,
                cells_fut,
                pure=False,
            )
            futures.append(fut)

        # one reduction task depending on all tile futures
        result_fut = client.submit(
            _reduce_parts, futures, len(cells), int(image_da.shape[0]), mode, scale, pure=False
        )
        return result_fut  # Future

    # ---- fallback: local/delayed path (OK for threaded scheduler) ----
    delayed_parts = []
    for (iy, ix, (xmin, xmax, ymin, ymax), idxs_) in tiles:
        bi = [0] * darr.ndim
        bi[c_ax] = 0
        bi[y_ax] = iy
        bi[x_ax] = ix
        chunk_ref = blocks[tuple(bi)]
        # local path can pass small lists directly; same semantics
        cells_subset = [cells[k] for k in idxs_]
        means_src = None
        if mode == "variance":
            if means_delayed is None:
                raise ValueError("means required for variance computation")
            means_src = means_delayed if not is_dask_collection(means_delayed) else dask.delayed(lambda x: x)(means_delayed)
        delayed_parts.append(
            dask.delayed(_aggregate_chunk_task)(
                chunk_ref, (xmin, xmax, ymin, ymax), idxs_, cells_subset, mode, means_src
            )
        )

    return dask.delayed(_reduce_parts)(delayed_parts, len(cells), int(image_da.shape[0]), mode, scale)


def to_intrinsic_from_image_coords(
    gdf_micron: gpd.GeoDataFrame,
    image_da,  # xarray.DataArray with coords "x" and "y" (monotonic, axis-aligned)
) -> gpd.GeoDataFrame:
    """
    Convert geometries from the image's world coords (micron) into the image's
    intrinsic pixel index space, without requiring the shapes to be in sdata.

    Works for axis-aligned grids with constant spacing (dx, dy), regardless of
    whether coords are ascending or descending.
    """
    x = image_da.coords["x"].values
    y = image_da.coords["y"].values

    # pixel spacing (allow descending coords)
    dx = float(np.mean(np.diff(x)))
    dy = float(np.mean(np.diff(y)))
    if dx == 0.0 or dy == 0.0:
        raise ValueError("Non-invertible coords: dx or dy is zero.")
    x0 = float(x[0])
    y0 = float(y[0])

    def world_to_pixel(geom):
        # shift world origin to (x0,y0), then scale to pixel units
        g = translate(geom, xoff=-x0, yoff=-y0)
        # scale by 1/dx, 1/dy; dy may be negative (common for images); that’s OK—
        # it flips the axis so 0,0 maps to the image’s first row/col.
        return scale(g, xfact=1.0/dx, yfact=1.0/dy, origin=(0, 0))

    out = gdf_micron.copy()
    out["geometry"] = [world_to_pixel(g) for g in out.geometry]
    return out

def aggregate_channels_lazy(
    sdata: SpatialData,
    image_key: str,
    shapes_gdf_micron: gpd.GeoDataFrame,
    expand_radius_ratio: float,
    mode: str,
    means_delayed=None,
):
    """
    Use your existing sopa utils (get_spatial_image, to_intrinsic, expand_radius),
    then call the parallel-safe aligned aggregator in lazy mode.
    """
    image = get_spatial_image(sdata, image_key)

    # convert micron -> intrinsic pixel coords against this image
    gdf_intrinsic = to_intrinsic_from_image_coords(shapes_gdf_micron, image)
    if expand_radius_ratio:
        gdf_intrinsic = expand_radius(gdf_intrinsic, expand_radius_ratio)

    return _aggregate_channels_aligned_parallel_lazy(
        image, gdf_intrinsic, mode, means_delayed=means_delayed
    )

def run_pipeline_parallel(
    args,                      # has data_path, recompute, etc.
    results_path: str | Path,
    compute_ficture_methods: list[str],
    logger,
    scheduler: str | None = None,
):
    """
    Full flow:
      1) prepare_ficture once
      2) build sdata images, set transforms
      3) build boundaries in parallel (micron space)
      4) batch compute area/mean across methods
      5) compute variance across methods (depends on each mean)
      6) write CSVs
    """
    results_path = Path(results_path) / "results"

    # 1) prepare_ficture ONCE
    logger.info("Preparing Ficture")
    stats = prepare_ficture(args.data_path, str(results_path), top_n_factors=1, logger=logger)
    area_covered = (stats["images"] > 0) #TODO: add seond prepare ficture
    images = stats["images"]
    del stats

    # 2) sdata images + transforms
    sdata = SpatialData()
    sdata["ficture_image_1"] = Image2DModel.parse(area_covered)
    sdata["ficture_image_2"] = Image2DModel.parse(images)
    del area_covered, images

    transformation = np.loadtxt(Path(args.data_path) / "images/micron_to_mosaic_pixel_transform.csv")
    aff_inv = Affine(transformation, input_axes=("x","y"), output_axes=("x","y")).inverse()
    set_transformation(sdata["ficture_image_1"], aff_inv, to_coordinate_system="global")
    set_transformation(sdata["ficture_image_2"], aff_inv, to_coordinate_system="global")

    # 3) boundaries in parallel (micron)
    logger.info("Building boundaries (parallel)")
    boundaries = build_boundaries_parallel(
        compute_ficture_methods, results_path, transformation, scheduler=scheduler
    )

    # 4) WAVE 1 — area (sum) + mean (average)
    logger.info("Wave 1: wiring area + mean")
    area_futs_or_delayed = []
    mean_futs_or_delayed = []
    keys = []

    for m, gdf_micron in boundaries.items():
        # (optional but recommended) rechunk both images so tiles are ~1024^2
        sdata["ficture_image_1"].data = sdata["ficture_image_1"].data.chunk({"c": -1, "y": 1024, "x": 1024})
        sdata["ficture_image_2"].data = sdata["ficture_image_2"].data.chunk({"c": -1, "y": 1024, "x": 1024})

        area_task = aggregate_channels_lazy(sdata, "ficture_image_1", gdf_micron, 0.0, "sum")
        mean_task = aggregate_channels_lazy(sdata, "ficture_image_2", gdf_micron, 0.0, "average")
        area_futs_or_delayed.append(area_task)
        mean_futs_or_delayed.append(mean_task)
        keys.append(m)

    if client is not None:
        progress(area_futs_or_delayed + mean_futs_or_delayed, notebook=True)
        area_vals = client.gather(area_futs_or_delayed)
        mean_vals = client.gather(mean_futs_or_delayed)
    else:
        from dask.diagnostics import ProgressBar
        with ProgressBar():
            area_vals = dask.compute(*area_futs_or_delayed, scheduler=scheduler)
            mean_vals = dask.compute(*mean_futs_or_delayed, scheduler=scheduler)

    areas = dict(zip(keys, area_vals))
    means = dict(zip(keys, mean_vals))

    # 5) WAVE 2 — variance (uses mean per method)
    logger.info("Wave 2: wiring variance")
    var_futs_or_delayed = []
    keys2 = []

    if client is not None:
        # scatter each method's mean once (small; ~MBs), submit variance futures
        for m, gdf_micron in boundaries.items():
            mean_ref = client.scatter(means[m], broadcast=True)
            task = aggregate_channels_lazy(
                sdata, "ficture_image_2", gdf_micron, 0.0, "variance", means_delayed=mean_ref
            )
            var_futs_or_delayed.append(task)
            keys2.append(m)
        progress(var_futs_or_delayed, notebook=True)
        var_vals = client.gather(var_futs_or_delayed)
    else:
        for m, gdf_micron in boundaries.items():
            task = aggregate_channels_lazy(
                sdata, "ficture_image_2", gdf_micron, 0.0, "variance", means_delayed=means[m]
            )
            var_futs_or_delayed.append(task)
            keys2.append(m)
        from dask.diagnostics import ProgressBar
        with ProgressBar():
            var_vals = dask.compute(*var_futs_or_delayed, scheduler=scheduler)

    vars_ = dict(zip(keys2, var_vals))

    # 6) write CSVs
    logger.info("Writing CSVs")
    for m, gdf_micron in boundaries.items():
        outdir = results_path / m / "Ficture_stats"
        outdir.mkdir(parents=True, exist_ok=True)
        index = gdf_micron.index.to_list()

        pd.DataFrame(
            means[m], index=index,
            columns=[f"fictureF21_{i}_mean_intensity_weighted" for i in range(21)]
        ).to_csv(outdir / "means_weight.csv")

        pd.DataFrame(
            vars_[m], index=index,
            columns=[f"fictureF21_{i}_variance_intensity_weighted" for i in range(21)]
        ).to_csv(outdir / "vars_weight.csv")

        pd.DataFrame(
            areas[m], index=index,
            columns=[f"fictureF21_{i}_area" for i in range(21)]
        ).to_csv(outdir / "area.csv")

    logger.info("Done.")