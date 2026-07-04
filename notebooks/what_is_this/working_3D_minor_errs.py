import logging
import os.path
from collections import Counter, defaultdict

import dask.bag as db
import geopandas as gpd
import numpy as np
import pyvista as pv
import shapely as shp
import vtkmodules.vtkInteractionStyle  # noqa
import vtkmodules.vtkRenderingOpenGL2  # noqa
from dask.distributed import Client
from scipy.spatial import cKDTree
from vtkmodules.vtkFiltersCore import vtkAppendPolyData
from vtkmodules.vtkFiltersModeling import vtkRuledSurfaceFilter

logger = logging.getLogger("shape_mapping")
logger.setLevel(logging.WARNING)
handler = logging.StreamHandler()  # or FileHandler(...)
handler.setFormatter(
    logging.Formatter("%(asctime)s [%(levelname)s] Cell %(cell_id)s: %(message)s")
)
logger.addHandler(handler)


def preparation_dome(seq, z_offsets, scale_factors, offset=0):
    xy = np.array(seq)[:, :2]
    z = np.unique(np.array(seq)[:, 2])
    assert len(z) == 1, f"{z}"
    n_pts = len(xy)
    center = xy.mean(axis=0)

    layers_mesh = []
    for z_off, scale in zip(z_offsets, scale_factors):
        scaled_xy = (xy - center) * scale + center
        layer = np.column_stack([scaled_xy, np.full(n_pts, z + z_off)])
        layers_mesh.append(layer)
    vertices = np.vstack(layers_mesh)
    faces = []

    for i in range(len(layers_mesh) - 1):
        lower = np.arange(i * n_pts, (i + 1) * n_pts)
        upper = np.arange((i + 1) * n_pts, (i + 2) * n_pts)

        for j in range(n_pts):
            a = lower[j]
            b = lower[(j + 1) % n_pts]
            c = upper[(j + 1) % n_pts]
            d = upper[j]
            faces.extend([4, a + offset, b + offset, c + offset, d + offset])
    return layers_mesh, vertices, faces


def determine_nearest_shape(point, trees):
    """Unchanged"""
    distances = [tree.query(point[:2])[0] for tree in trees]
    return int(np.argmin(distances))


def determine_shapes_mapping(
    bottom_list, top_list, *, cell_id=None, bottom_layer=None, top_layer=None
):
    """Returns (relationships, conflicts).
    Logs each conflict at WARNING, including cell_id, bottom_layer and top_layer in the message.
    """
    # Build KD-trees
    bottom_trees = [cKDTree(pts[:, :2]) for pts in bottom_list]
    top_trees = [cKDTree(pts[:, :2]) for pts in top_list]

    # Helper: per-vertex vote with centroid fallback
    def vote_with_tiebreak(coords, trees):
        votes = [determine_nearest_shape(pt, trees) for pt in coords]
        cnt = Counter(votes)
        top2 = cnt.most_common(2)
        if len(top2) == 1 or top2[0][1] != top2[1][1]:
            return top2[0][0]
        xy = coords[:, :2]
        poly = shp.Polygon(xy)
        cx, cy = poly.centroid.x, poly.centroid.y
        # feed only 2D point into the KD-tree
        return determine_nearest_shape((cx, cy), trees)

    # 1) bottom→top decisions
    decisions_bottom = [vote_with_tiebreak(shape, top_trees) for shape in bottom_list]
    # 2) top→bottom decisions
    decisions_top = [vote_with_tiebreak(shape, bottom_trees) for shape in top_list]

    # Invert to dicts
    bottom_to_top = {i: t for i, t in enumerate(decisions_bottom)}
    top_to_bottoms = defaultdict(list)
    for b, t in bottom_to_top.items():
        top_to_bottoms[t].append(b)

    bottom_to_tops = defaultdict(list)
    for t, b in enumerate(decisions_top):
        bottom_to_tops[b].append(t)

    # Classify
    relationships = []
    conflicts = []
    used_b, used_t = set(), set()

    # Pass 1: many-to-one (but detect many-to-many)
    for t, bottoms in top_to_bottoms.items():
        if len(bottoms) <= 1:
            continue
        multi = [b for b in bottoms if len(bottom_to_tops[b]) > 1]
        if multi:
            for b in multi:
                entry = {
                    "type": "many-to-many",
                    "bottom": b,
                    "top": t,
                    "bottom_layer": bottom_layer,
                    "top_layer": top_layer,
                }
                conflicts.append(entry)
                logger.warning(
                    f"[layers {bottom_layer}→{top_layer}] "
                    f"many-to-many conflict bottom={b} ↔ top={t}",
                    extra={"cell_id": cell_id},
                )
        else:
            relationships.append(
                {
                    "type": "many-to-one",
                    "bottom": bottoms,
                    "top": t,
                    "bottom_layer": bottom_layer,
                    "top_layer": top_layer,
                }
            )
            used_t.add(t)
            used_b.update(bottoms)

    # Pass 2: one-to-many
    for b, tops in bottom_to_tops.items():
        if b in used_b or len(tops) <= 1:
            continue
        if any(t in used_t for t in tops):
            for t in tops:
                entry = {
                    "type": "many-to-many",
                    "bottom": b,
                    "top": t,
                    "bottom_layer": bottom_layer,
                    "top_layer": top_layer,
                }
                conflicts.append(entry)
                logger.warning(
                    f"[layers {bottom_layer}→{top_layer}] "
                    f"many-to-many conflict bottom={b} ↔ top={t}",
                    extra={"cell_id": cell_id},
                )
        else:
            relationships.append(
                {
                    "type": "one-to-many",
                    "bottom": b,
                    "top": tops,
                    "bottom_layer": bottom_layer,
                    "top_layer": top_layer,
                }
            )
            used_b.add(b)
            used_t.update(tops)

    # Pass 3: one-to-one
    for b, t in bottom_to_top.items():
        if b in used_b or t in used_t:
            continue
        if len(top_to_bottoms[t]) == 1 and len(bottom_to_tops[b]) == 1:
            relationships.append(
                {
                    "type": "one-to-one",
                    "bottom": b,
                    "top": t,
                    "bottom_layer": bottom_layer,
                    "top_layer": top_layer,
                }
            )
            used_b.add(b)
            used_t.add(t)

    # Leftovers → unclassified
    for b, t in bottom_to_top.items():
        if b not in used_b or t not in used_t:
            entry = {
                "type": "unclassified",
                "bottom": b,
                "top": t,
                "bottom_layer": bottom_layer,
                "top_layer": top_layer,
            }
            conflicts.append(entry)
            logger.warning(
                f"[layers {bottom_layer}→{top_layer}] "
                f"unclassified mapping bottom={b} → top={t}",
                extra={"cell_id": cell_id},
            )

    return relationships, conflicts


def extract_exterior_rings(row):
    """Given a GeoPandas row (pd.Series) with .geometry (Polygon/MultiPolygon)
    and .ZLevel, returns a list of (n_pts×3) numpy arrays for each exterior ring.
    """
    geom = row.Geometry
    z = row.ZLevel
    rings = []

    parts = [geom] if isinstance(geom, shp.Polygon) else list(geom.geoms)
    for poly in parts:
        coords2d = np.array(poly.exterior.coords)  # shape (n_pts, 2)
        zs = np.full((coords2d.shape[0], 1), z)  # shape (n_pts, 1)
        rings.append(np.hstack([coords2d, zs]))  # shape (n_pts, 3)

    return rings


def ruled_surface_between(bottom_pts, top_pts, n_samples=None):
    """Build a ruled surface between two closed xyz loops by:
    1) creating two PolyData polylines,
    2) appending them into one dataset,
    3) feeding that to vtkRuledSurfaceFilter.
    """

    # 1) make PyVista polylines
    def make_polyline(pts):
        # 🔥 CHANGED: explicitly close the loop by repeating first point
        if not np.allclose(pts[0], pts[-1]):
            pts = np.vstack([pts, pts[0]])
        n = pts.shape[0]
        pd = pv.PolyData(pts)
        # now lines=[n, 0,1,…,(n-1)] will connect last→first as well
        pd.lines = np.hstack([[n], np.arange(n)])
        return pd

    pd_bot = make_polyline(bottom_pts)
    pd_top = make_polyline(top_pts)

    # 2) append into one PolyData
    app = vtkAppendPolyData()
    app.AddInputData(pd_bot)  # both lines go into port 0
    app.AddInputData(pd_top)
    app.Update()
    combined = pv.wrap(app.GetOutput())

    # 3) run the filter on the combined dataset
    rsf = vtkRuledSurfaceFilter()
    rsf.SetInputData(combined)  # 🔥 CHANGED: only one input
    if n_samples is None:
        n_samples = max(len(bottom_pts), len(top_pts))
    rsf.SetResolution(n_samples, 1)
    rsf.SetRuledModeToResample()
    rsf.Update()

    return pv.wrap(rsf.GetOutput()).triangulate()


def align_loops(bottom_pts, top_pts):
    """Rotate & optionally reverse top_pts so that:
      - top_pts[0] is the closest to bottom_pts[0]
      - both loops have the same orientation (sign of signed area).
    Returns a new array top_aligned of same shape.
    """
    # 1) find nearest top index to bottom_pts[0]
    tree = cKDTree(top_pts[:, :2])
    _, idx0 = tree.query(bottom_pts[0, :2])
    top = np.roll(top_pts, -idx0, axis=0)

    # 2) ensure same winding (signed area in XY)
    def signed_area(pts):
        x, y = pts[:, 0], pts[:, 1]
        return 0.5 * np.dot(x, np.roll(y, -1) - np.roll(y, 1))

    if np.sign(signed_area(bottom_pts)) != np.sign(signed_area(top)):
        top = top[::-1]

    return top


def triangulate(cell_geometry, z_offsets, scale_factors, cell_id):
    faces = []
    vertices = np.empty((0, 3))

    # --- BUILD BOTTOM CAP & BOTTOM RINGS ---
    bottom_row = cell_geometry[cell_geometry.ZLevel == cell_geometry.ZLevel.min()].iloc[
        0
    ]
    bottom_rings = extract_exterior_rings(bottom_row)

    # 1a) bottom dome caps
    for ring in bottom_rings:
        layers, dome_verts, dome_faces = preparation_dome(
            ring, -z_offsets, scale_factors
        )
        faces.extend(dome_faces)
        vertices = np.vstack([vertices, dome_verts])

    # 1) collect all rings per layer
    z_vals = sorted(cell_geometry.ZLevel.unique())
    layers = []
    for z in z_vals:
        row = cell_geometry[cell_geometry.ZLevel == z].iloc[0]
        layers.append(extract_exterior_rings(row))

    ring_bases = []  # list of lists: one sublist per layer
    for rings in layers:
        layer_bases = []
        for ring in rings:
            layer_bases.append(vertices.shape[0])
            vertices = np.vstack([vertices, ring])
        ring_bases.append(layer_bases)

    # … later, inside your loop over layer‐pairs …
    for idx in range(len(layers) - 1):
        bot_rings = layers[idx]
        top_rings = layers[idx + 1]
        bot_z, top_z = z_vals[idx], z_vals[idx + 1]

        rels, confs = determine_shapes_mapping(
            bot_rings,
            top_rings,
            cell_id=cell_id,
            bottom_layer=bot_z,
            top_layer=top_z,
        )
        for rel in rels:
            b_idxs = (
                rel["bottom"] if isinstance(rel["bottom"], list) else [rel["bottom"]]
            )
            t_idxs = rel["top"] if isinstance(rel["top"], list) else [rel["top"]]
            for bi in b_idxs:
                for ti in t_idxs:
                    # 1) grab the two full rings
                    bot_pts = bot_rings[bi]
                    top_pts = top_rings[ti]

                    # 🔥 CHANGED: if many→one or one→many, downsample the big ring
                    if rel["type"] in ("many-to-one", "one-to-many"):
                        # decide which is smaller
                        if len(bot_pts) <= len(top_pts):
                            small, big = bot_pts, top_pts
                            small_is_bot = True
                        else:
                            small, big = top_pts, bot_pts
                            small_is_bot = False

                        N = len(small)
                        # pick N nearest from the big ring
                        tree = cKDTree(big[:, :2])
                        _, idxs = tree.query(small[:, :2])
                        sampled_big = big[idxs]

                        # close both loops by appending the first point
                        def close(pts):
                            return (
                                np.vstack([pts, pts[0]])
                                if not np.allclose(pts[0], pts[-1])
                                else pts
                            )

                        small_closed = close(small)
                        big_closed = close(sampled_big)

                        # reassign in the correct order
                        if small_is_bot:
                            bot_pts, top_pts = small_closed, big_closed
                        else:
                            bot_pts, top_pts = big_closed, small_closed

                    # 🔥 CHANGED: align start‐points & winding
                    top_pts = align_loops(bot_pts, top_pts)

                    # 2) build the ruled surface
                    wall = ruled_surface_between(bot_pts, top_pts)

                    # 3) append its vertices & faces into your global arrays
                    base = vertices.shape[0]
                    pts = wall.points
                    tris = wall.faces.reshape(-1, 4)[:, 1:]

                    vertices = np.vstack([vertices, pts])
                    for tri in tris:
                        faces.extend([3, *(tri + base)])

    # --- FINALLY TOP CAP (similar to bottom) ---
    # --- 4) Top dome cap (🔥 CHANGED: correctly placed AFTER all side‐walls) ---
    top_row = cell_geometry[cell_geometry.ZLevel == z_vals[-1]].iloc[0]
    top_rings = extract_exterior_rings(top_row)
    for ring in top_rings:
        offset = vertices.shape[0]  # 🔥 CHANGED
        layers, dome_verts, dome_faces = preparation_dome(
            ring,
            z_offsets,
            scale_factors,
            offset=offset,  # 🔥 CHANGED
        )
        faces.extend(dome_faces)  # 🔥 CHANGED
        vertices = np.vstack([vertices, dome_verts])

    mesh = pv.PolyData(vertices, np.array(faces))
    mesh.clean(inplace=True)
    hole_size = 1e5  # adjust this to slightly larger than your largest hole radius
    mesh = mesh.fill_holes(hole_size)
    return mesh


def build_3D_cell(cell_geometry, id, z_cap_height=1.4, z_samples=10):
    z_offsets = np.linspace(0, z_cap_height, z_samples)
    scale_factors = np.sqrt(np.clip(1 - (z_offsets / z_cap_height) ** 2, 0, 1))
    mesh = triangulate(cell_geometry, z_offsets, scale_factors, id)
    path = os.path.join(save_path, "tmp", f"cell_{id}.vtp")
    mesh.save(path)
    return {"id": id, "results": path}


def aggregate_dict(dicts):
    return [x for x in dicts]


save_path = "/Users/jonasflor/Desktop/test_Cellpose3D_Paraview"


def aggregate_list(partitions):
    block = pv.MultiBlock()
    for dicts in partitions:
        for x in dicts:
            block.append(pv.load(x["results"]), name=str(x["id"]))
    block.save(save_path)
    return save_path


def main():
    layers = gpd.read_parquet(
        "/Users/jonasflor/Downloads/cellpose2_mosaic_space.parquet"
    )
    layers = layers[[not x.is_empty for x in layers.geometry]]

    def shapely_to_pyvista(shapely_geometry):
        # Check if it's a multipolygon
        if isinstance(shapely_geometry, shp.geometry.MultiPolygon):
            coords = []
            for poly in list(shapely_geometry.geoms):
                coords_new = poly.exterior.coords
                coords.append(coords_new)
            # Return a list of meshes for the multipolygon
            return coords
        else:
            # If it's just a single polygon
            coords = shapely_geometry.exterior.coords
            return coords

    layers["coords"] = [shapely_to_pyvista(x) for x in layers.geometry]
    layers["ZLevel"] = layers.ZLevel * 10

    cell_groups = []
    for cell_id, subdf in layers.groupby("EntityID"):
        cell_groups.append((subdf, cell_id))
    pv.OFF_SCREEN = True
    client = Client(threads_per_worker=1, n_workers=10)
    print(client.dashboard_link)
    b = db.from_sequence(cell_groups, npartitions=40)
    delayed = b.map(
        lambda cid: build_3D_cell(cid[0], cid[1], z_cap_height=14)
    ).reduction(aggregate_dict, aggregate_list)
    delayed.compute()


if __name__ == "__main__":
    main()
