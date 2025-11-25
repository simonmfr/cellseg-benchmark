from pathlib import Path

import numpy as np
from geopandas import GeoDataFrame
from scipy.ndimage import binary_fill_holes, convolve
from shapely.geometry import GeometryCollection, MultiPolygon, Point, Polygon, box
from shapely.ops import unary_union
from shapely.prepared import prep
from shapely.strtree import STRtree
from skimage.measure import find_contours
from skimage.measure import label as cc_label

try:
    import shapely
    from packaging.version import Version
    from shapely.prepared import prep as _prepare

    _SHAPELY2 = Version(shapely.__version__).major >= 2
except Exception:
    from shapely.prepare import prepare as _prepare

    _SHAPELY2 = False

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon as MplPolygon


# ------------------------- 0) Build Grid -----------------------------------
def gridify(adata, banksy_clustering):
    # --- 0) Pull data once ---
    sp = np.asarray(adata.obsm["spatial_microns"], dtype=float)  # columns: [x, y] in µm
    labels = np.asarray(adata.obs[banksy_clustering])

    # --- 1) Infer geo with TOP-LEFT physical origin ---
    minx, maxx = sp[:, 0].min(), sp[:, 0].max()
    miny, maxy = sp[:, 1].min(), sp[:, 1].max()

    ux = np.unique(sp[:, 0])
    uy = np.unique(sp[:, 1])
    dx = np.median(np.diff(ux)) if ux.size > 1 else 1.0
    dy = np.median(np.diff(uy)) if uy.size > 1 else 1.0

    x0 = float(minx)
    y0 = float(maxy)  # TOP edge in physical coords
    geo = (x0, y0, float(dx), float(dy))

    # --- 2) Grid extent (include both ends) ---
    ext_x = int(np.rint((maxx - x0) / dx)) + 1
    ext_y = int(np.rint((y0 - miny) / dy)) + 1  # note y0 - miny (top-to-bottom)

    # --- 3) Vectorized binning (x,y) -> (row,col) for UPPER-LEFT origin ---
    cols = np.rint((sp[:, 0] - x0) / dx).astype(int)  # 0..ext_x-1
    rows = np.rint((y0 - sp[:, 1]) / dy).astype(int)  # 0..ext_y-1 (top down)

    # --- 4) In-bounds mask ---
    inb = (cols >= 0) & (cols < ext_x) & (rows >= 0) & (rows < ext_y)
    rows, cols, labels = rows[inb], cols[inb], labels[inb]

    # --- 5) Enforce uniqueness (one point per pixel) ---
    key = rows.astype(np.int64) * ext_x + cols.astype(np.int64)
    uk, cnt = np.unique(key, return_counts=True)
    if np.any(cnt > 1):
        dup_rc = np.column_stack(np.unravel_index(uk[cnt > 1], (ext_y, ext_x)))
        raise ValueError(
            f"Multiple cells map to the same pixel. First duplicate (row,col): {dup_rc[0].tolist()} "
            f"(use an aggregation if this is expected)."
        )

    # --- 6) Build grid ---
    grid = -np.ones((ext_y, ext_x), dtype=np.int32)
    grid[rows, cols] = labels
    return grid, geo


# ---------- 1) CLEAN SMALL HOLES (ENCLOSED INCLUSIONS) ---------------------


def relabel_small_holes(grid, *, min_hole_area_um2, dx, dy, connectivity=1):
    """grid: 2D array of label *codes* (-1 for background is fine, but not required)
    Relabels fully enclosed inclusions smaller than threshold to the surrounding majority label.
    """
    G = np.asarray(grid).copy()
    pix_area = float(dx) * float(dy)
    min_hole_px = int(np.floor(min_hole_area_um2 / pix_area + 1e-9))

    K8 = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], int)  # for boundary majority

    for L in np.unique(G[G >= 0]):
        maskL = G == L

        filled = binary_fill_holes(maskL)  # fills all holes of L
        holes = filled & (~maskL)
        if not holes.any():
            continue

        lab = cc_label(holes, connectivity=connectivity)
        for hid in range(1, lab.max() + 1):
            H = lab == hid
            if H.sum() > min_hole_px:
                continue  # keep larger internal structure

            boundary = convolve(H.astype(int), K8, mode="constant", cval=0) > 0
            neigh = G[boundary & (~H)]
            neigh = neigh[neigh >= 0]
            if neigh.size == 0:
                continue

            vals, cnts = np.unique(neigh, return_counts=True)
            maj = int(vals[np.argmax(cnts)])
            if maj == int(L):
                G[H] = L
    return G


# ---------- 2) (OPTIONAL) DROP TINY ISLANDS OF EACH LABEL ------------------


def remove_small_islands(grid, *, min_island_area_um2, dx, dy, connectivity=1):
    """Reassigns small connected components of each label to the local boundary majority."""
    G = np.asarray(grid).copy()
    pix_area = float(dx) * float(dy)
    min_comp_px = int(np.floor(min_island_area_um2 / pix_area + 1e-9))

    K8 = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], int)
    for L in np.unique(G[G >= 0]):
        lab = cc_label(G == L, connectivity=connectivity)
        for cid in range(1, lab.max() + 1):
            comp = lab == cid
            sz = int(comp.sum())
            if sz == 0 or sz > min_comp_px:
                continue
            boundary = convolve(comp.astype(int), K8, mode="constant", cval=0) > 0
            neigh = G[boundary & (~comp)]
            neigh = neigh[neigh >= 0]
            if neigh.size == 0:
                G[comp] = -1
                continue
            vals, cnts = np.unique(neigh, return_counts=True)
            newL = int(vals[np.argmax(cnts)])
            if newL != int(L):
                G[comp] = newL
    return G


# ---------- 3) POLYGONIZE EACH CONNECTED COMPONENT (NO HULLS) --------------


def polygons_per_component_from_grid(clean_grid, geo, value_map=None, connectivity=1):
    """clean_grid: 2D array of label *codes* (>=0); -1 treated as background.
    geo: (x0, y0, dx, dy) with TOP-LEFT origin; pixel (row,col) -> (x0+col*dx, y0 - row*dy)
    """
    x0, y0, dx, dy = [float(v) for v in geo]
    regions = []
    codes = np.unique(clean_grid[clean_grid >= 0])

    for code in codes:
        comp_lab = cc_label(clean_grid == code, connectivity=connectivity)
        for cid in range(1, comp_lab.max() + 1):
            reg = (comp_lab == cid).astype(float)
            curves = find_contours(reg, level=0.5)
            if not curves:
                continue

            rings = []
            for crv in curves:
                rr, cc = crv[:, 0], crv[:, 1]
                x = x0 + cc * dx
                y = y0 - rr * dy  # <— top-left origin mapping
                a = np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))  # signed area
                rings.append((a, np.column_stack([x, y])))

            rings.sort(key=lambda t: abs(t[0]), reverse=True)
            exterior = rings[0][1]
            holes = [r[1] for r in rings[1:] if r[0] < 0]

            poly = Polygon(exterior, holes if holes else None)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if poly.is_empty:
                continue

            regions.append(
                {
                    "poly": poly,
                    "code": int(code),
                    "value": (
                        value_map[int(code)] if value_map is not None else int(code)
                    ),
                    "comp_id": int(cid - 1),
                }
            )
    return regions


# ---------- 4) SHAPELY INDEX + POINT → REGION MAPPING ----------------------


def build_region_index(regions):
    geoms = [r["poly"] for r in regions]
    tree = STRtree(geoms)
    prepared = [prep(g) for g in geoms]
    if _SHAPELY2:

        def candidates(geom):
            return tree.query(geom)  # idx array
    else:
        id_to_idx = {id(g): i for i, g in enumerate(geoms)}

        def candidates(geom):
            return [id_to_idx[id(g)] for g in tree.query(geom)]

    return prepared, candidates


def map_points_to_regions(
    points_xy, regions, prepared, candidates_fn, *, include_boundary=True
):
    pts = np.asarray(points_xy, float)
    region_idx = np.full(len(pts), -1, int)
    cluster_val = [None] * len(pts)
    comp_id = np.full(len(pts), -1, int)

    for i, (x, y) in enumerate(pts):
        p = Point(float(x), float(y))
        for k in candidates_fn(p):
            ok = prepared[k].covers(p) if include_boundary else prepared[k].contains(p)
            if ok:
                region_idx[i] = k
                cluster_val[i] = regions[k]["value"]
                comp_id[i] = regions[k]["comp_id"]
                break
    return region_idx, cluster_val, comp_id


# ---------- 5) PIPELINE ENTRY POINT (FROM YOUR EXISTING GRID) --------------


def segment_from_grid_and_map_points(
    label_grid,  # 2D array of label codes (>=0) and/or background (-1)
    geo,  # (x0, y0, dx, dy); here dx=dy=25 for your case
    points_xy,  # (N,2) points in same coord frame as geo
    *,
    min_hole_area_um2,  # e.g. 2000.0
    min_island_area_um2=None,  # optional, e.g. 2000.0
    connectivity=1,  # 1 (4-neigh) preserves thin gaps; 2 (8-neigh) smoother
    value_map=None,  # optional mapping code->original label value
):
    # 1) hole-fix
    x0, y0, dx, dy = geo
    clean = relabel_small_holes(
        label_grid,
        min_hole_area_um2=min_hole_area_um2,
        dx=dx,
        dy=dy,
        connectivity=connectivity,
    )

    # 2) optional tiny-island cleanup
    if min_island_area_um2 is not None and min_island_area_um2 > 0:
        clean = remove_small_islands(
            clean,
            min_island_area_um2=min_island_area_um2,
            dx=dx,
            dy=dy,
            connectivity=connectivity,
        )

    # 3) polygons per connected component
    regions = polygons_per_component_from_grid(
        clean, geo, value_map=value_map, connectivity=connectivity
    )

    # 4) spatial index + mapping
    prepared, cand = build_region_index(regions)
    mapping = map_points_to_regions(
        points_xy, regions, prepared, cand, include_boundary=True
    )

    return clean, regions, mapping


def polygons_per_component_exact(clean_grid, geo, value_map=None, connectivity=1):
    x0, y0, dx, dy = [float(v) for v in geo]
    regions = []
    for code in np.unique(clean_grid[clean_grid >= 0]):
        comp = cc_label(clean_grid == code, connectivity=connectivity)
        for cid in range(1, comp.max() + 1):
            rr, cc = np.where(comp == cid)
            if rr.size == 0:
                continue
            # TOP-LEFT origin: row r spans [y_top, y_bot] = [y0 - r*dy, y0 - (r+1)*dy]
            polys = []
            for r, c in zip(rr, cc):
                x_left = x0 + c * dx
                x_right = x0 + (c + 1) * dx
                y_top = y0 - r * dy
                y_bot = y0 - (r + 1) * dy
                ymin, ymax = (y_bot, y_top) if y_bot < y_top else (y_top, y_bot)
                polys.append(box(x_left, ymin, x_right, ymax))

            poly = unary_union(polys)
            if not poly.is_empty:
                regions.append(
                    {
                        "poly": poly,
                        "code": int(code),
                        "value": (
                            value_map[int(code)] if value_map is not None else int(code)
                        ),
                        "comp_id": int(cid - 1),
                    }
                )
    return regions


# ---------- plot_steps.py ----------
# ======== small helpers ========


def _extent_from_geo(grid, geo):
    """Imshow extent in physical units from (x0, y0, dx, dy) with TOP-LEFT origin."""
    H, W = grid.shape
    x0, y0, dx, dy = [float(v) for v in geo]
    # left, right, bottom, top — note bottom < top
    return [x0, x0 + W * dx, y0 - H * dy, y0]


def _show_grid(grid, geo, title="", ax=None, cmap="tab20", vmin=None, vmax=None):
    """Generic label-grid renderer (physical coordinates; TOP-LEFT origin)."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.figure
    extent = _extent_from_geo(grid, geo)
    im = ax.imshow(
        grid,
        origin="upper",
        extent=extent,
        cmap=cmap,
        interpolation="nearest",
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(title)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    plt.colorbar(im, ax=ax)
    return fig, ax


def _resolve_highlight_indices(regions, highlight):
    """Return a list of region indices to highlight based on `highlight` spec."""
    if highlight is None:
        return []
    # numeric index
    if isinstance(highlight, int):
        return [highlight] if 0 <= highlight < len(regions) else []

    # tuple (code, comp_id) or (code, None)
    if isinstance(highlight, tuple) and len(highlight) == 2:
        code, comp = highlight
        idxs = [
            i
            for i, r in enumerate(regions)
            if r["code"] == code and (comp is None or r["comp_id"] == comp)
        ]
        return idxs

    # dict form
    if isinstance(highlight, dict):
        if "index" in highlight:
            i = int(highlight["index"])
            return [i] if 0 <= i < len(regions) else []
        code = highlight.get("code", None)
        comp = highlight.get("comp_id", None)
        if code is not None:
            return [
                i
                for i, r in enumerate(regions)
                if r["code"] == code and (comp is None or r["comp_id"] == comp)
            ]

    # unsupported → nothing
    return []


def _overlay_polygons(
    regions,
    ax=None,
    *,
    highlight=None,
    base_edge_kw=None,
    highlight_edge_kw=None,
    highlight_face_alpha=0.20,
):
    """Draw polygon boundaries for a list of region dicts.
    Optionally highlight one or more regions.

    Parameters
    ----------
    regions : list of dicts
        Output of polygons_per_component_*; each has keys: 'poly','code','value','comp_id'.
    ax : matplotlib Axes, optional
    highlight : int | (code, comp_id) | (code, None) | dict
        Region selector. See _resolve_highlight_indices() doc above.
    base_edge_kw : dict, optional
        Matplotlib kwargs for non-highlighted boundaries (e.g., {'linewidth':1.2, 'alpha':0.9}).
    highlight_edge_kw : dict, optional
        Matplotlib kwargs for highlighted boundary (e.g., {'linewidth':2.5}).
    highlight_face_alpha : float
        If >0, fills the highlighted region(s) with the edge color at this alpha.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))

    # defaults
    if base_edge_kw is None:
        base_edge_kw = {"linewidth": 1.5, "alpha": 1.0}
    if highlight_edge_kw is None:
        highlight_edge_kw = {"linewidth": 3.0}

    # 1) draw ALL boundaries (thin)
    lines = []
    for r in regions:
        g = r["poly"]
        geoms = [g] if g.geom_type == "Polygon" else list(g.geoms)
        for poly in geoms:
            x, y = poly.exterior.xy
            lines.append(np.column_stack([x, y]))
            for hole in poly.interiors:
                hx, hy = hole.xy
                lines.append(np.column_stack([hx, hy]))
    if lines:
        lc = LineCollection(lines, **base_edge_kw)
        ax.add_collection(lc)

    # 2) overlay HIGHLIGHT region(s) (thick + optional fill)
    hi_idxs = _resolve_highlight_indices(regions, highlight)
    if hi_idxs:
        for i in hi_idxs:
            g = regions[i]["poly"]
            geoms = [g] if g.geom_type == "Polygon" else list(g.geoms)
            for poly in geoms:
                # edge
                x, y = poly.exterior.xy
                ax.plot(x, y, **highlight_edge_kw)
                # holes (draw edges)
                for hole in poly.interiors:
                    hx, hy = hole.xy
                    ax.plot(hx, hy, **highlight_edge_kw)
                # face (optional, semi-transparent)
                if highlight_face_alpha and "color" in highlight_edge_kw:
                    face_color = highlight_edge_kw["color"]
                else:
                    face_color = None
                if face_color and highlight_face_alpha > 0:
                    patch = MplPolygon(
                        np.column_stack([x, y]),
                        closed=True,
                        facecolor=face_color,
                        alpha=float(highlight_face_alpha),
                        edgecolor="none",
                    )
                    ax.add_patch(patch)

    return ax


# ======== step A: original ========


def plot_original_grid(label_grid, geo, *, cmap="tab20"):
    """Input:
      label_grid: (H,W) int array (codes; -1 allowed)
      geo: (x0, y0, dx, dy) in µm
    Returns: fig, ax
    """
    fig, ax = _show_grid(label_grid, geo, title="A) Original grid", cmap=cmap)
    return fig, ax


# ======== step B: after hole relabeling ========


def plot_after_holes(
    label_grid, geo, *, relabel_func, min_hole_area_um2, connectivity=1, cmap="tab20"
):
    """Runs the provided relabel_func and plots the result.
    Inputs:
      relabel_func: callable(grid, min_hole_area_um2, dx, dy, connectivity) -> grid
    Returns: clean_holes_grid, fig, ax
    """
    x0, y0, dx, dy = geo
    clean_holes = relabel_func(
        label_grid,
        min_hole_area_um2=min_hole_area_um2,
        dx=dx,
        dy=dy,
        connectivity=connectivity,
    )
    fig, ax = _show_grid(
        clean_holes,
        geo,
        title=f"B) After hole relabeling (≤ {min_hole_area_um2:g} µm²)",
        cmap=cmap,
    )
    return clean_holes, fig, ax


# ======== step C: after tiny-island removal (optional) ========


def plot_after_islands(
    clean_holes_grid,
    geo,
    *,
    remove_islands_func,
    min_island_area_um2,
    connectivity=1,
    cmap="tab20",
):
    """Runs the provided remove_islands_func and plots the result.
    Inputs:
      remove_islands_func: callable(grid, min_island_area_um2, dx, dy, connectivity) -> grid
    Returns: clean_final_grid, fig, ax
    """
    x0, y0, dx, dy = geo
    clean_final = remove_islands_func(
        clean_holes_grid,
        min_island_area_um2=min_island_area_um2,
        dx=dx,
        dy=dy,
        connectivity=connectivity,
    )
    fig, ax = _show_grid(
        clean_final,
        geo,
        title=f"C) After tiny-island removal (≤ {min_island_area_um2:g} µm²)",
        cmap=cmap,
    )
    return clean_final, fig, ax


# ======== step D: final boundaries (+ optional points) ========


def plot_final_boundaries(
    clean_grid,
    geo,
    *,
    polygonize_func,
    connectivity=1,
    value_map=None,
    points_xy=None,
    points_kwargs=None,
    cmap="tab20",
):
    """Builds polygons from the cleaned grid and overlays them.
    Inputs:
      polygonize_func: callable(clean_grid, geo, value_map, connectivity) -> list of region dicts
      points_xy: optional (N,2) array in µm to overlay
      points_kwargs: dict of matplotlib scatter kwargs (size, marker, etc.)
    Returns: regions, fig, ax
    """
    regions = polygonize_func(
        clean_grid, geo, value_map=value_map, connectivity=connectivity
    )
    fig, ax = _show_grid(clean_grid, geo, title="D) Final boundaries", cmap=cmap)
    _overlay_polygons(regions, ax)
    if points_xy is not None:
        kw = {"s": 20, "marker": "x", "linewidths": 1.0}
        if points_kwargs:
            kw.update(points_kwargs)
        pts = np.asarray(points_xy, float)
        ax.scatter(pts[:, 0], pts[:, 1], **kw)
    return regions, fig, ax


def plot_progress_pipeline(
    adata,
    sample,
    banksy_clustering,
    *,
    relabel_func,  # function(grid, min_hole_area_um2, dx, dy, connectivity) -> grid
    remove_islands_func=None,  # function(grid, min_island_area_um2, dx, dy, connectivity) -> grid
    polygonize_func=None,  # function(clean_grid, geo, value_map, connectivity) -> regions
    points_xy=None,  # optional (N,2) in µm; plotted on final panel
    min_hole_area_um2=2500.0,
    min_island_area_um2=None,
    connectivity=1,
    value_map=None,
    figsize=(14, 9),
    cmap="tab20",
    save=None,  # None => don't save; str/Path => save figure there and don't return fig
    save_dpi=200,
    save_kwargs=None,
):
    """End-to-end visualizer:
      Panel A: original grid
      Panel B: after hole relabeling
      Panel C: after tiny-island removal (if requested)
      Panel D: final boundaries + (optional) points

    Returns:
      clean_holes, clean_final, regions
    """
    adata_tmp = adata[adata.obs["sample"] == sample].copy()
    label_grid, geo = gridify(adata_tmp, banksy_clustering)

    x0, y0, dx, dy = geo
    # A) original
    fig, axs = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)
    _show_grid(label_grid, geo, title="A) Original grid", ax=axs[0, 0], cmap=cmap)

    # B) after holes
    clean_holes = relabel_func(
        label_grid,
        min_hole_area_um2=min_hole_area_um2,
        dx=dx,
        dy=dy,
        connectivity=connectivity,
    )
    _show_grid(
        clean_holes, geo, title="B) After hole relabeling", ax=axs[0, 1], cmap=cmap
    )

    # C) after islands (optional)
    if remove_islands_func is not None and (
        min_island_area_um2 is not None and min_island_area_um2 > 0
    ):
        clean_final = remove_islands_func(
            clean_holes,
            min_island_area_um2=min_island_area_um2,
            dx=dx,
            dy=dy,
            connectivity=connectivity,
        )
        c_title = f"C) After tiny-island removal (≤ {min_island_area_um2:.0f} µm²)"
    else:
        clean_final = clean_holes
        c_title = "C) (skipped) Tiny-island removal"
    _show_grid(clean_final, geo, title=c_title, ax=axs[1, 0], cmap=cmap)

    # D) boundaries + points
    axD = axs[1, 1]
    _show_grid(
        clean_final, geo, title="D) Final boundaries + points", ax=axD, cmap=cmap
    )
    regions = []
    if polygonize_func is not None:
        regions = polygonize_func(
            clean_final, geo, value_map=value_map, connectivity=connectivity
        )
        _overlay_polygons(regions, ax=axD)
    if points_xy is not None:
        _overlay_points(points_xy, ax=axD, size=18, marker="x")

    if save is not None:
        p = Path(save)
        p.parent.mkdir(parents=True, exist_ok=True)
        kw = {"bbox_inches": "tight", "dpi": save_dpi}
        if save_kwargs:
            kw.update(save_kwargs)
        fig.savefig(p, **kw)
        plt.close(fig)
        return clean_holes, clean_final, regions

    return clean_holes, clean_final, regions, fig, geo


def _overlay_points(points_xy, ax=None, size=12, marker="o"):
    """Scatter arbitrary points (µm) on current axes."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))
    pts = np.asarray(points_xy, float)
    ax.scatter(pts[:, 0], pts[:, 1], s=size, marker=marker)
    return ax


def _flatten_to_polys(geom):
    """Return a list of Polygon(s) from any Shapely geometry, ignoring non-areas."""
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return list(geom.geoms)
    if isinstance(geom, GeometryCollection):
        out = []
        for g in geom.geoms:
            if isinstance(g, (Polygon, MultiPolygon)):
                out.extend(_flatten_to_polys(g))
        return out
    return []  # lines/points/etc. are ignored


def dissolve_slide_labels(
    label_polys_for_slide: dict,
    *,
    join_distance: float = 0.0,  # µm; 0 = only edge-sharing/overlap
    min_area: float = 0.0,  # µm²; drop tiny slivers after dissolve
) -> dict:
    """Input:  { label: [Polygon, Polygon, ...], ... }
    Output: { label: [MergedPolygon, MergedPolygon, ...] }  (each item is a connected piece)
    - If join_distance > 0: close small gaps by buffering out + union + buffer back.
    - min_area removes tiny artefacts produced by buffering or topology fixes.
    """
    merged = {}
    for lab, polys in label_polys_for_slide.items():
        polys = [p for p in polys if (p is not None and not p.is_empty)]
        if not polys:
            merged[lab] = []
            continue

        if join_distance > 0:
            # bridge small gaps, then shrink back
            grown = [p.buffer(join_distance) for p in polys]
            u = unary_union(grown).buffer(-join_distance)
        else:
            u = unary_union(polys)

        # repair minor invalidities
        try:
            if not u.is_valid:
                u = u.buffer(0)
        except Exception:
            pass

        parts = [
            g for g in _flatten_to_polys(u) if (min_area <= 0 or g.area >= min_area)
        ]
        merged[lab] = parts
    return merged


# ---------------- utilities ----------------
def _bounds_from_polys(label_polys):
    # label_polys: dict[label -> list[Polygon]]
    xmin = ymin = np.inf
    xmax = ymax = -np.inf
    for polys in label_polys.values():
        for g in polys:
            if g.is_empty:
                continue
            bxmin, bymin, bxmax, bymax = g.bounds
            xmin = min(xmin, bxmin)
            ymin = min(ymin, bymin)
            xmax = max(xmax, bxmax)
            ymax = max(ymax, bymax)
    if not np.isfinite([xmin, ymin, xmax, ymax]).all():
        xmin = ymin = 0.0
        xmax = ymax = 1.0
    return xmin, ymin, xmax, ymax


def _label_color_map(labels, cmap_name="tab20"):
    # deterministic color per label string
    uniq = sorted(set(labels))
    cmap = cm.get_cmap(cmap_name, max(10, len(uniq)))
    return {lab: cmap(i % cmap.N) for i, lab in enumerate(uniq)}


def _iter_all_polylines(geom):
    # yield exterior line then hole lines as Nx2 arrays
    if geom.geom_type == "Polygon":
        geoms = [geom]
    else:
        geoms = list(geom.geoms)
    for poly in geoms:
        x, y = poly.exterior.xy
        yield np.column_stack([x, y])
        for hole in poly.interiors:
            hx, hy = hole.xy
            yield np.column_stack([hx, hy])


# -------------- main plotting ----------------


def plot_slide_regions(
    regions_by_slide: dict,
    slide_id: str,
    *,
    # optional background
    label_grid=None,  # 2D int array; if given, shown as background
    geo=None,  # (x0,y0,dx,dy) needed if label_grid is given
    cmap_grid="tab20",
    show_grid_colorbar: bool = False,
    # styling
    base_edge_kw=None,  # kwargs for non-highlighted outlines
    highlight=None,  # label string OR tuple ("label", index) to highlight
    highlight_edge_kw=None,  # kwargs for highlighted outline (e.g. {"color":"red","linewidth":3})
    highlight_face_alpha=0.18,  # translucent fill for highlighted region(s)
    dissolve: bool = True,
    join_distance: float = 0.0,  # µm, how far to “bridge” gaps within a label
    min_area: float = 0.0,
    # figure output
    figsize=(8, 7),
    save=None,  # None -> return fig,ax; str/Path -> save and return None
    save_dpi=300,
    save_kwargs=None,
):
    """Plot borders for a single slide from your nested dict:
      regions_by_slide[slide_id] = { label: [Polygon, Polygon, ...], ... }

    highlight:
      - pass a label string: e.g. "CTX" -> highlight all CTX polygons
      - pass ("CTX", 2)     -> highlight only the 3rd CTX polygon
    """
    if slide_id not in regions_by_slide:
        raise KeyError(f"slide_id '{slide_id}' not in regions_by_slide")

    label_polys: dict[str, list] = regions_by_slide[slide_id]

    # pick colors per label
    label_colors = _label_color_map(label_polys.keys())

    # defaults
    if base_edge_kw is None:
        base_edge_kw = {"linewidth": 3.4, "alpha": 1.0, "color": (0, 0, 0, 0.75)}
    if highlight_edge_kw is None:
        highlight_edge_kw = {"linewidth": 6.0, "color": "crimson"}

    if slide_id not in regions_by_slide:
        raise KeyError(f"slide_id '{slide_id}' not in regions_by_slide")
    label_polys_src: dict[str, list] = regions_by_slide[slide_id]

    # 2) (DO THIS NOW) Optionally dissolve per label
    label_polys = (
        dissolve_slide_labels(
            label_polys_src, join_distance=join_distance, min_area=min_area
        )
        if dissolve
        else label_polys_src
    )

    # 3) Set up figure/axes
    fig, ax = plt.subplots(figsize=figsize)

    # 4) Background (optional)
    if label_grid is not None and geo is not None:
        extent = _extent_from_geo(label_grid, geo)
        im = ax.imshow(
            label_grid,
            origin="upper",
            extent=extent,  # <-- was "lower"
            cmap=cmap_grid,
            interpolation="nearest",
        )
        if show_grid_colorbar:
            plt.colorbar(im, ax=ax)
    else:
        # Use (dissolved) polygons to set bounds
        xmin, ymin, xmax, ymax = _bounds_from_polys(label_polys)
        pad_x = 0.02 * (xmax - xmin) if xmax > xmin else 1.0
        pad_y = 0.02 * (ymax - ymin) if ymax > ymin else 1.0
        ax.set_xlim(xmin - pad_x, xmax + pad_x)
        ax.set_ylim(ymin - pad_y, ymax + pad_y)

    # draw all boundaries first
    for lab, polys in label_polys.items():
        if not polys:
            continue
        col = label_colors.get(lab, (0, 0, 0, 1))
        lines = []
        for g in polys:
            if g.is_empty:
                continue
            for seg in _iter_all_polylines(g):
                lines.append(seg)
        if lines:
            lc = LineCollection(lines, **base_edge_kw)
            lc.set_color(col)
            ax.add_collection(lc)

    # highlight logic
    if highlight is not None:
        if isinstance(highlight, tuple) and len(highlight) == 2:
            lab_sel, idx = highlight
            polys = label_polys.get(lab_sel, [])
            to_highlight = [
                (lab_sel, i, p)
                for i, p in enumerate(polys)
                if (idx is None or i == idx)
            ]
        else:
            # treat as label name -> highlight all in that label
            lab_sel = str(highlight)
            polys = label_polys.get(lab_sel, [])
            to_highlight = [(lab_sel, i, p) for i, p in enumerate(polys)]

        for lab_sel, i, g in to_highlight:
            if g is None or g.is_empty:
                continue
            # thick edge
            for seg in _iter_all_polylines(g):
                ax.plot(seg[:, 0], seg[:, 1], **highlight_edge_kw)
            # translucent fill for the exterior only
            try:
                geoms = [g] if g.geom_type == "Polygon" else list(g.geoms)
                face_color = highlight_edge_kw.get("color", None)
                if face_color and highlight_face_alpha > 0:
                    for poly in geoms:
                        x, y = poly.exterior.xy
                        patch = MplPolygon(
                            np.column_stack([x, y]),
                            closed=True,
                            facecolor=face_color,
                            alpha=float(highlight_face_alpha),
                            edgecolor="none",
                        )
                        ax.add_patch(patch)
            except Exception:
                pass

    # legend (labels only if not too many)
    if len(label_polys) <= 20:
        handles = []
        for lab, col in label_colors.items():
            h = plt.Line2D([0], [0], color=col, lw=2, label=lab)
            handles.append(h)
        ax.legend(
            handles=handles,
            bbox_to_anchor=(1.02, 1.0),
            loc="upper left",
            borderaxespad=0.0,
        )

    ax.set_aspect("equal")
    ax.set_title(f"Borders — {slide_id}")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")

    # save or return fig
    if save is not None:
        p = Path(save)
        p.parent.mkdir(parents=True, exist_ok=True)
        kw = {"dpi": save_dpi, "bbox_inches": "tight"}
        if save_kwargs:
            kw.update(save_kwargs)
        fig.savefig(p, **kw)
        plt.close(fig)
        return

    return fig, ax


# ---------- build per-slide spatial index ----------
def _build_slide_index(label_polys_for_slide: dict):
    """Returns:
    tree        : STRtree over all polygons of this slide
    prepared    : list of prepared geometries (same order as tree input)
    meta        : list of tuples (label, poly_idx_within_label)
    cand_fn     : function(Point)->iter(indices) (handles Shapely 1.8 vs 2.x)
    """
    geoms, meta = [], []
    for lab, polys in label_polys_for_slide.items():
        for j, g in enumerate(polys):
            if g is None or g.is_empty:
                continue
            geoms.append(g)
            meta.append((lab, j))

    if not geoms:
        # empty index
        dummy = [Polygon()]
        tree = STRtree(dummy)
        prepared = [_prepare(dummy[0])]

        def cand_fn(_):
            return []

        return tree, prepared, meta, cand_fn

    tree = STRtree(geoms)
    prepared = [_prepare(g) for g in geoms]
    if _SHAPELY2:

        def cand_fn(geom):  # returns integer indices
            return tree.query(geom)
    else:
        # Shapely 1.8 returns geometries; map to indices
        id_to_idx = {id(g): i for i, g in enumerate(geoms)}

        def cand_fn(geom):
            return [id_to_idx[id(g)] for g in tree.query(geom)]

    return tree, prepared, meta, cand_fn


# ---------- map points for ONE slide ----------
def _map_points_for_slide(
    points_xy: np.ndarray,
    tree,
    prepared,
    meta,
    cand_fn,
    *,
    include_boundary: bool = True,
):
    """Returns:
    labels_out : length-N object array with label (or None)
    poly_idx   : length-N int array with polygon index within that label (or -1)
    """
    pts = np.asarray(points_xy, float)
    N = len(pts)
    labels_out = np.empty(N, dtype=object)
    labels_out[:] = None
    poly_idx = np.full(N, -1, dtype=int)

    for i, (x, y) in enumerate(pts):
        p = Point(float(x), float(y))
        for k in cand_fn(p):
            ok = prepared[k].covers(p) if include_boundary else prepared[k].contains(p)
            if ok:
                lab, j = meta[k]
                labels_out[i] = lab
                poly_idx[i] = j
                break
    return labels_out, poly_idx


# ---------- AnnData-aware wrapper ----------
def map_points_to_regions_from_anndata(
    adata,
    regions_by_slide: dict,
    *,
    coord_key: str = "spatial_microns",  # obsm key with (N,2) microns
    slide_key: str = "sample",  # obs key with slide id per cell
    slides: list | None = None,  # restrict to subset; None = all present in adata
    include_boundary: bool = True,  # covers vs contains
    dissolve: bool = False,  # merge touching polys per label first
    join_distance: float = 0.0,  # µm; only used if dissolve=True
    min_area: float = 0.0,  # µm²; only used if dissolve=True
    index_kind: str = "name",  # 'name' -> use obs_names, 'pos' -> integer positions
    return_df: bool = True,
):
    """For each slide in `adata.obs[slide_key]`, map points in `adata.obsm[coord_key]`
    to region labels using polygons in `regions_by_slide[slide_id][label] = [Polygon,...]`.

    Returns dict[slide_id] with:
      - "indices_and_labels": list[(obs_id, label_or_None)]
      - "labels": np.ndarray (N,)
      - "poly_index": np.ndarray (N,)
      - if return_df: "df": pandas.DataFrame with columns [obs_id, slide, x, y, label, poly_index]
    """
    # --- sanity
    if coord_key not in adata.obsm_keys():
        raise KeyError(f"obsm['{coord_key}'] not found.")
    if slide_key not in adata.obs:
        raise KeyError(f"obs['{slide_key}'] not found.")

    coords_all = np.asarray(adata.obsm[coord_key], float)
    slide_series = adata.obs[slide_key].astype(str)

    # which slides?
    if slides is None:
        slides = sorted(slide_series.unique())
    else:
        slides = [str(s) for s in slides]

    results = {}
    for sid in slides:
        mask = (slide_series == sid).values
        if not mask.any():
            continue
        if sid not in regions_by_slide:
            # No polygons for this slide → all None
            idxs = np.nonzero(mask)[0]
            labels_out = np.array([None] * mask.sum(), dtype=object)
            poly_idx = np.full(mask.sum(), -1, dtype=int)
            obs_ids = adata.obs_names.values[idxs] if index_kind == "name" else idxs
            idx_and_lab = list(zip(obs_ids, labels_out.tolist()))
            out = {
                "indices_and_labels": idx_and_lab,
                "labels": labels_out,
                "poly_index": poly_idx,
            }
            if return_df:
                import pandas as pd

                pts = coords_all[mask, :]
                out["df"] = pd.DataFrame(
                    {
                        "obs_id": obs_ids,
                        "slide": sid,
                        "x": pts[:, 0],
                        "y": pts[:, 1],
                        "label": labels_out,
                        "poly_index": poly_idx,
                    }
                )
            results[sid] = out
            continue

        # gather points for this slide
        pts = coords_all[mask, :]
        idxs = np.nonzero(mask)[0]
        obs_ids = adata.obs_names.values[idxs] if index_kind == "name" else idxs

        # dissolve (optional) then index
        label_polys_src = regions_by_slide[sid]
        label_polys = (
            dissolve_slide_labels(
                label_polys_src, join_distance=join_distance, min_area=min_area
            )
            if dissolve
            else label_polys_src
        )

        tree, prepared, meta, cand_fn = _build_slide_index(label_polys)
        labels_out, poly_idx = _map_points_for_slide(
            pts, tree, prepared, meta, cand_fn, include_boundary=include_boundary
        )

        idx_and_lab = list(zip(obs_ids, labels_out.tolist()))
        out = {
            "indices_and_labels": idx_and_lab,
            "labels": labels_out,
            "poly_index": poly_idx,
        }

        if return_df:
            import pandas as pd

            out["df"] = pd.DataFrame(
                {
                    "obs_id": obs_ids,
                    "slide": sid,
                    "x": pts[:, 0],
                    "y": pts[:, 1],
                    "label": labels_out,
                    "poly_index": poly_idx,
                }
            )

        results[sid] = out

    return results


def to_gdf(d):
    rows = []
    for sample, labels in d.items():
        for label, geoms in labels.items():
            for i, g in enumerate(geoms):
                rows.append({"sample": sample, "label": label, "idx": i, "geometry": g})
    return GeoDataFrame(rows, geometry="geometry")
