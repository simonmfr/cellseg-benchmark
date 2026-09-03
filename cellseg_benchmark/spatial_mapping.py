"""Turn a raster spatial-domain clustering into region polygons and map cells into them.

The pipeline is: label raster -> morphological cleanup -> one polygon per
connected component -> point-in-polygon assignment of cells.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from geopandas import GeoDataFrame
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon as MplPolygon
from scipy.ndimage import binary_fill_holes, convolve
from shapely.geometry import Point, box
from shapely.ops import unary_union
from shapely.prepared import prep
from shapely.strtree import STRtree
from skimage.measure import label as cc_label

# 8-connected structuring element, used to read the labels bordering a component.
_K8 = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], int)


# --------------------------------------------------------------------------- #
# 0) raster
# --------------------------------------------------------------------------- #
def gridify(adata, label_key, coord_key="spatial_microns"):
    """Bin cells of one sample onto the regular raster they came from.

    Args:
        adata: AnnData of a single sample, with integer labels in `.obs`.
        label_key: `.obs` column holding integer label codes (>= 0).
        coord_key: `.obsm` key with (N, 2) coordinates in microns.

    Returns:
        Tuple of the (H, W) int32 label grid (-1 = background) and
        `geo = (x0, y0, dx, dy)` with a TOP-LEFT physical origin.
    """
    sp = np.asarray(adata.obsm[coord_key], dtype=float)
    labels = np.asarray(adata.obs[label_key])

    minx, maxx = sp[:, 0].min(), sp[:, 0].max()
    miny, maxy = sp[:, 1].min(), sp[:, 1].max()

    ux, uy = np.unique(sp[:, 0]), np.unique(sp[:, 1])
    dx = np.median(np.diff(ux)) if ux.size > 1 else 1.0
    dy = np.median(np.diff(uy)) if uy.size > 1 else 1.0

    x0, y0 = float(minx), float(maxy)  # top-left corner in physical coords
    geo = (x0, y0, float(dx), float(dy))

    ext_x = int(np.rint((maxx - x0) / dx)) + 1
    ext_y = int(np.rint((y0 - miny) / dy)) + 1

    cols = np.rint((sp[:, 0] - x0) / dx).astype(int)
    rows = np.rint((y0 - sp[:, 1]) / dy).astype(int)

    inb = (cols >= 0) & (cols < ext_x) & (rows >= 0) & (rows < ext_y)
    rows, cols, labels = rows[inb], cols[inb], labels[inb]

    key = rows.astype(np.int64) * ext_x + cols.astype(np.int64)
    uk, cnt = np.unique(key, return_counts=True)
    if np.any(cnt > 1):
        dup = np.column_stack(np.unravel_index(uk[cnt > 1], (ext_y, ext_x)))
        raise ValueError(
            f"Multiple cells map to the same pixel, first at (row, col) "
            f"{dup[0].tolist()}. gridify expects a raster segmentation."
        )

    grid = -np.ones((ext_y, ext_x), dtype=np.int32)
    grid[rows, cols] = labels
    return grid, geo


# --------------------------------------------------------------------------- #
# 1) cleanup
# --------------------------------------------------------------------------- #
def relabel_small_holes(grid, *, min_hole_area_um2, dx, dy, connectivity=1):
    """Relabel fully enclosed inclusions smaller than the threshold.

    Args:
        grid: 2D array of label codes (-1 = background).
        min_hole_area_um2: Holes at or below this area are filled.
        dx: Pixel width in microns.
        dy: Pixel height in microns.
        connectivity: 1 = 4-neighbour, 2 = 8-neighbour.

    Returns:
        A copy of `grid` with small enclosed holes assigned to the surrounding label.
    """
    g = np.asarray(grid).copy()
    min_hole_px = int(np.floor(min_hole_area_um2 / (float(dx) * float(dy)) + 1e-9))

    for lab in np.unique(g[g >= 0]):
        mask = g == lab
        holes = binary_fill_holes(mask) & (~mask)
        if not holes.any():
            continue

        cc = cc_label(holes, connectivity=connectivity)
        for hid in range(1, cc.max() + 1):
            hole = cc == hid
            if hole.sum() > min_hole_px:
                continue  # keep genuine internal structure
            border = convolve(hole.astype(int), _K8, mode="constant", cval=0) > 0
            neigh = g[border & (~hole)]
            neigh = neigh[neigh >= 0]
            if neigh.size == 0:
                continue
            vals, cnts = np.unique(neigh, return_counts=True)
            if int(vals[np.argmax(cnts)]) == int(lab):
                g[hole] = lab
    return g


def remove_small_islands(grid, *, min_island_area_um2, dx, dy, connectivity=1):
    """Reassign small connected components of each label to their boundary majority.

    Args:
        grid: 2D array of label codes (-1 = background).
        min_island_area_um2: Components at or below this area are dissolved.
        dx: Pixel width in microns.
        dy: Pixel height in microns.
        connectivity: 1 = 4-neighbour, 2 = 8-neighbour.

    Returns:
        A copy of `grid` without small islands.
    """
    g = np.asarray(grid).copy()
    min_comp_px = int(np.floor(min_island_area_um2 / (float(dx) * float(dy)) + 1e-9))

    for lab in np.unique(g[g >= 0]):
        cc = cc_label(g == lab, connectivity=connectivity)
        for cid in range(1, cc.max() + 1):
            comp = cc == cid
            if comp.sum() == 0 or comp.sum() > min_comp_px:
                continue
            border = convolve(comp.astype(int), _K8, mode="constant", cval=0) > 0
            neigh = g[border & (~comp)]
            neigh = neigh[neigh >= 0]
            if neigh.size == 0:
                g[comp] = -1
                continue
            vals, cnts = np.unique(neigh, return_counts=True)
            new = int(vals[np.argmax(cnts)])
            if new != int(lab):
                g[comp] = new
    return g


# --------------------------------------------------------------------------- #
# 2) polygonize
# --------------------------------------------------------------------------- #
def polygons_per_component_exact(clean_grid, geo, value_map=None, connectivity=1):
    """Build one polygon per connected component, following exact pixel edges.

    Args:
        clean_grid: 2D array of label codes (-1 = background).
        geo: `(x0, y0, dx, dy)` with TOP-LEFT origin.
        value_map: Optional mapping from code to the value stored per region.
        connectivity: 1 = 4-neighbour, 2 = 8-neighbour.

    Returns:
        List of dicts with keys `poly`, `code`, `value` and `comp_id`.
    """
    x0, y0, dx, dy = (float(v) for v in geo)
    regions = []
    for code in np.unique(clean_grid[clean_grid >= 0]):
        cc = cc_label(clean_grid == code, connectivity=connectivity)
        for cid in range(1, cc.max() + 1):
            rr, ccols = np.where(cc == cid)
            if rr.size == 0:
                continue
            polys = [
                box(x0 + c * dx, y0 - (r + 1) * dy, x0 + (c + 1) * dx, y0 - r * dy)
                for r, c in zip(rr, ccols)
            ]
            poly = unary_union(polys)
            if poly.is_empty:
                continue
            regions.append(
                {
                    "poly": poly,
                    "code": int(code),
                    "value": value_map[int(code)] if value_map else int(code),
                    "comp_id": int(cid - 1),
                }
            )
    return regions


# --------------------------------------------------------------------------- #
# 3) map cells into regions
# --------------------------------------------------------------------------- #
def _build_slide_index(label_polys):
    """Build an STRtree over all polygons of one slide.

    Args:
        label_polys: Mapping from label to a list of polygons.

    Returns:
        Tuple of prepared geometries, `(label, index_within_label)` metadata and a
        callable returning candidate indices for a query geometry.
    """
    geoms, meta = [], []
    for lab, polys in label_polys.items():
        for j, g in enumerate(polys):
            if g is not None and not g.is_empty:
                geoms.append(g)
                meta.append((lab, j))

    if not geoms:
        return [], meta, lambda _: []

    tree = STRtree(geoms)
    return [prep(g) for g in geoms], meta, tree.query


def _map_points_for_slide(points_xy, prepared, meta, cand_fn, *, include_boundary=True):
    """Assign each point of one slide to a region.

    Args:
        points_xy: (N, 2) coordinates in microns.
        prepared: Prepared geometries from `_build_slide_index`.
        meta: Metadata from `_build_slide_index`.
        cand_fn: Candidate lookup from `_build_slide_index`.
        include_boundary: Use `covers` instead of `contains`.

    Returns:
        Tuple of an object array of labels (None where unassigned) and an int
        array of polygon indices within the label (-1 where unassigned).
    """
    pts = np.asarray(points_xy, float)
    labels = np.empty(len(pts), dtype=object)
    labels[:] = None
    poly_idx = np.full(len(pts), -1, dtype=int)

    for i, (x, y) in enumerate(pts):
        p = Point(float(x), float(y))
        for k in cand_fn(p):
            hit = prepared[k].covers(p) if include_boundary else prepared[k].contains(p)
            if hit:
                labels[i], poly_idx[i] = meta[k]
                break
    return labels, poly_idx


def map_points_to_regions_from_anndata(
    adata,
    regions_by_slide: dict,
    *,
    coord_key: str = "spatial_microns",
    slide_key: str = "sample",
    slides: list | None = None,
    include_boundary: bool = True,
    index_kind: str = "name",
    return_df: bool = True,
):
    """Map every cell of every slide onto the region polygons of its slide.

    Args:
        adata: AnnData with coordinates in `.obsm[coord_key]`.
        regions_by_slide: `{slide: {label: [Polygon, ...]}}`.
        coord_key: `.obsm` key with (N, 2) coordinates in microns.
        slide_key: `.obs` column holding the slide id.
        slides: Restrict to these slides. None uses all slides present.
        include_boundary: Assign points lying exactly on a border.
        index_kind: `"name"` for `obs_names`, `"pos"` for integer positions.
        return_df: Also return a tidy DataFrame per slide.

    Returns:
        Dict keyed by slide with `labels`, `poly_index`, `indices_and_labels` and,
        if requested, `df`.
    """
    if coord_key not in adata.obsm:
        raise KeyError(f"obsm['{coord_key}'] not found.")
    if slide_key not in adata.obs:
        raise KeyError(f"obs['{slide_key}'] not found.")

    coords_all = np.asarray(adata.obsm[coord_key], float)
    slide_series = adata.obs[slide_key].astype(str)
    slides = (
        sorted(slide_series.unique()) if slides is None else [str(s) for s in slides]
    )

    results = {}
    for sid in slides:
        mask = (slide_series == sid).values
        if not mask.any():
            continue

        idxs = np.nonzero(mask)[0]
        obs_ids = adata.obs_names.values[idxs] if index_kind == "name" else idxs
        pts = coords_all[mask, :]

        if sid in regions_by_slide:
            prepared, meta, cand_fn = _build_slide_index(regions_by_slide[sid])
            labels, poly_idx = _map_points_for_slide(
                pts, prepared, meta, cand_fn, include_boundary=include_boundary
            )
        else:
            labels = np.array([None] * mask.sum(), dtype=object)
            poly_idx = np.full(mask.sum(), -1, dtype=int)

        out = {
            "indices_and_labels": list(zip(obs_ids, labels.tolist())),
            "labels": labels,
            "poly_index": poly_idx,
        }
        if return_df:
            out["df"] = pd.DataFrame(
                {
                    "obs_id": obs_ids,
                    "slide": sid,
                    "x": pts[:, 0],
                    "y": pts[:, 1],
                    "label": labels,
                    "poly_index": poly_idx,
                }
            )
        results[sid] = out
    return results


def to_gdf(regions_by_slide: dict) -> GeoDataFrame:
    """Flatten `{slide: {label: [Polygon, ...]}}` into a GeoDataFrame.

    Args:
        regions_by_slide: Nested region dict.

    Returns:
        GeoDataFrame with columns `sample`, `label`, `idx` and `geometry`.
    """
    rows = [
        {"sample": sample, "label": label, "idx": i, "geometry": g}
        for sample, labels in regions_by_slide.items()
        for label, geoms in labels.items()
        for i, g in enumerate(geoms)
    ]
    return GeoDataFrame(rows, geometry="geometry")


# --------------------------------------------------------------------------- #
# 4) inspection plots
# --------------------------------------------------------------------------- #
def _extent_from_geo(grid, geo):
    """Return the imshow extent in microns for a TOP-LEFT origin grid.

    Args:
        grid: 2D label grid.
        geo: `(x0, y0, dx, dy)`.

    Returns:
        `[left, right, bottom, top]`.
    """
    h, w = grid.shape
    x0, y0, dx, dy = (float(v) for v in geo)
    return [x0, x0 + w * dx, y0 - h * dy, y0]


def _show_grid(grid, geo, title="", ax=None, cmap="tab20"):
    """Render a label grid in physical coordinates.

    Args:
        grid: 2D label grid.
        geo: `(x0, y0, dx, dy)`.
        title: Axis title.
        ax: Existing axis, or None to create one.
        cmap: Matplotlib colormap.

    Returns:
        Tuple of figure and axis.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(
        grid,
        origin="upper",
        extent=_extent_from_geo(grid, geo),
        cmap=cmap,
        interpolation="nearest",
    )
    ax.set_title(title)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    plt.colorbar(im, ax=ax)
    return ax.figure, ax


def _resolve_highlight_indices(regions, highlight):
    """Resolve a highlight selector to region indices.

    Args:
        regions: List of region dicts.
        highlight: An index, a `(code, comp_id)` tuple (comp_id may be None), or None.

    Returns:
        List of region indices.
    """
    if highlight is None:
        return []
    if isinstance(highlight, int):
        return [highlight] if 0 <= highlight < len(regions) else []
    if isinstance(highlight, tuple) and len(highlight) == 2:
        code, comp = highlight
        return [
            i
            for i, r in enumerate(regions)
            if r["code"] == code and (comp is None or r["comp_id"] == comp)
        ]
    return []


def _overlay_polygons(regions, ax=None, *, highlight=None, highlight_edge_kw=None):
    """Draw region boundaries, optionally highlighting a selection.

    Args:
        regions: List of region dicts.
        ax: Existing axis, or None to create one.
        highlight: Selector understood by `_resolve_highlight_indices`.
        highlight_edge_kw: Matplotlib kwargs for the highlighted outline.

    Returns:
        The axis.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))
    highlight_edge_kw = highlight_edge_kw or {"linewidth": 3.0, "color": "red"}

    lines = []
    for r in regions:
        g = r["poly"]
        for poly in [g] if g.geom_type == "Polygon" else g.geoms:
            lines.append(np.column_stack(poly.exterior.xy))
            lines.extend(np.column_stack(h.xy) for h in poly.interiors)
    if lines:
        ax.add_collection(LineCollection(lines, linewidth=1.5, color=(0, 0, 0, 0.75)))

    for i in _resolve_highlight_indices(regions, highlight):
        g = regions[i]["poly"]
        for poly in [g] if g.geom_type == "Polygon" else g.geoms:
            xy = np.column_stack(poly.exterior.xy)
            ax.plot(xy[:, 0], xy[:, 1], **highlight_edge_kw)
            ax.add_patch(
                MplPolygon(
                    xy,
                    closed=True,
                    facecolor=highlight_edge_kw.get("color", "red"),
                    alpha=0.2,
                    edgecolor="none",
                )
            )
    return ax


def plot_cleanup_pipeline(
    adata,
    sample,
    label_key,
    *,
    min_hole_area_um2,
    min_island_area_um2=None,
    connectivity=1,
    value_map=None,
    save=None,
    cmap="tab20",
):
    """Show the raster before and after cleanup, plus the resulting borders.

    Use this to tune `min_hole_area_um2` and `min_island_area_um2` per sample.

    Args:
        adata: AnnData containing all samples.
        sample: Sample to plot.
        label_key: `.obs` column with integer label codes.
        min_hole_area_um2: Hole-filling threshold.
        min_island_area_um2: Island-removal threshold, or None to skip.
        connectivity: 1 = 4-neighbour, 2 = 8-neighbour.
        value_map: Optional mapping from code to region value.
        save: Path to save the figure, or None to return it.
        cmap: Matplotlib colormap.

    Returns:
        Tuple of the cleaned grid, the regions and `geo`.
    """
    sub = adata[adata.obs["sample"] == sample]
    grid, geo = gridify(sub, label_key)
    _, _, dx, dy = geo

    fig, axs = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
    _show_grid(grid, geo, title=f"{sample}: raw", ax=axs[0], cmap=cmap)

    clean = relabel_small_holes(
        grid,
        min_hole_area_um2=min_hole_area_um2,
        dx=dx,
        dy=dy,
        connectivity=connectivity,
    )
    if min_island_area_um2:
        clean = remove_small_islands(
            clean,
            min_island_area_um2=min_island_area_um2,
            dx=dx,
            dy=dy,
            connectivity=connectivity,
        )
    _show_grid(clean, geo, title=f"{sample}: cleaned", ax=axs[1], cmap=cmap)

    regions = polygons_per_component_exact(
        clean, geo, value_map=value_map, connectivity=connectivity
    )
    _show_grid(clean, geo, title=f"{sample}: borders", ax=axs[2], cmap=cmap)
    _overlay_polygons(regions, ax=axs[2])

    if save is not None:
        p = Path(save)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=200, bbox_inches="tight")
        plt.close(fig)
    return clean, regions, geo
