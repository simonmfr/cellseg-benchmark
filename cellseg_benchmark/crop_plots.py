"""Segmentation crops: find candidate regions and plot cell boundaries of methods in a region."""

from functools import cache
from pathlib import Path

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spatialdata as sd
import spatialdata_plot  # noqa: F401
from matplotlib.colors import LinearSegmentedColormap, to_rgb

from . import _constants
from . import ficture_utils as fu
from .metrics.ovrl import low_z_integrity_per_cell

dapi_cmap = LinearSegmentedColormap.from_list("dapi", ["white", "#4747ff"])


def _adata_path(sample, method):
    return Path(_constants.BASE_PATH) / "analysis" / sample.split("_")[0] / method / "adatas/adata_integrated.h5ad.gz"


def find_candidates(cohort, exclude=(), tile=145, quantile=0.85, n=8, ref="vpt_2D_DAPI_PolyT"):
    """Tiles inside tissue near a density quantile, highest low z integrity share first; at most 2 per sample, 4 tiles apart."""
    with h5py.File(_adata_path(cohort, ref)) as f:
        cells = pd.DataFrame(f["obsm/spatial"][:], columns=["x", "y"],
                             index=np.asarray(ad.io.read_elem(f["obs"][f["obs"].attrs["_index"]])).astype(str)).assign(
            sample=np.asarray(ad.io.read_elem(f["obs/sample"])).astype(str))
    cells = cells[~cells["sample"].isin(exclude)]
    cells["low_z"] = low_z_integrity_per_cell(ad.AnnData(obs=cells[["sample"]]), ref, _constants.BASE_PATH).to_numpy()
    tiles = (cells.assign(x0=cells.x // tile * tile, y0=cells.y // tile * tile)
             .groupby(["sample", "x0", "y0"]).agg(n=("low_z", "size"), low_z=("low_z", "mean")).reset_index())
    full = set(tiles.loc[tiles.n >= tiles.n.median(), ["sample", "x0", "y0"]].itertuples(index=False, name=None))
    inside = [all((s, x + dx * tile, y + dy * tile) in full for dx in (-1, 0, 1) for dy in (-1, 0, 1))
              for s, x, y in tiles[["sample", "x0", "y0"]].itertuples(index=False)]
    picked = []
    for r in tiles[inside].assign(d=(tiles.n - tiles.n.quantile(quantile)).abs()).sort_values(["d", "low_z"], ascending=[True, False]).itertuples():
        same = [p for p in picked if p.sample == r.sample]
        if len(same) < 2 and all(np.hypot(r.x0 - p.x0, r.y0 - p.y0) >= 4 * tile for p in same):
            picked.append(r)
        if len(picked) == n:
            break
    return pd.DataFrame(picked).drop(columns="Index")


@cache
def _cell_types(sample, method):
    with h5py.File(_adata_path(sample, method)) as f:
        mask = ad.io.read_elem(f["obs/sample"]) == sample
        ids = pd.Index(ad.io.read_elem(f["obs"][f["obs"].attrs["_index"]]))[mask]
        return ids, np.asarray(ad.io.read_elem(f["obs/cell_type_revised"]))[mask].astype(str)


def draw(ax, sdata, sample, method, box, style="outlines", qc="post", lw=0.4):
    """Draw one method in box (xmin, ymin, size in µm) as outlines on DAPI, filled cell types or outlines on FICTURE transcripts.

    qc="both" adds QC-removed cells as grey outlines.
    """
    xlim, ylim = (box[0], box[0] + box[2]), (box[1], box[1] + box[2])
    crop = sd.SpatialData.init_from_elements({k: sd.bounding_box_query(
        sdata[k], axes=["x", "y"], min_coordinate=[xlim[0], ylim[0]], max_coordinate=[xlim[1], ylim[1]],
        target_coordinate_system="micron") for k in [f"{sample}_z3", f"boundaries_{method}"]})
    if style == "outlines":
        crop.pl.render_images(f"{sample}_z3", channel="DAPI", cmap=dapi_cmap, colorbar=False, scale="full").pl.show(
            coordinate_systems="micron", ax=ax, title="")
    if style == "transcripts":
        t = pd.read_parquet(Path(_constants.BASE_PATH) / "samples" / sample / "results/Ficture/ficture_transcript_factors.parquet",
                            columns=["x", "y", "factor"],
                            filters=[("x", ">=", xlim[0]), ("x", "<=", xlim[1]), ("y", ">=", ylim[0]), ("y", "<=", ylim[1])])
        ct = t.factor.map(fu.factor_cell_types(fu.find_ficture_output(sample, _constants.BASE_PATH)))
        ax.scatter(t.x, t.y, s=1.5 * lw, c=ct.map(_constants.cell_type_colors).fillna("#d9d9d9"), lw=0, alpha=0.6, rasterized=True)
    g = sd.transform(crop[f"boundaries_{method}"], to_coordinate_system="micron")
    if method.startswith(tuple(_constants.methods_3D)):
        g = g.reset_index(drop=True)[["cell_id", "geometry"]].dissolve("cell_id")
    gi = g.index.astype(str)
    ids, ct = _cell_types(sample, method)
    ct = pd.Series(ct, index=ids.where(ids.isin(gi), ids.str.rsplit("-", n=1).str[0]))
    ct = ct[~ct.index.duplicated()].reindex(gi).to_numpy()
    if style == "celltypes":
        fill = pd.Series(ct).map(_constants.cell_type_colors).fillna("#e0e0e0").to_numpy()
        order = np.argsort(-g.area.to_numpy())
        g.iloc[order].plot(ax=ax, aspect=None, color=fill[order].tolist(), lw=0.4, edgecolor=[tuple(0.7 * v for v in to_rgb(c)) for c in fill[order]])
    else:
        kept = pd.notna(ct)
        if qc == "both":
            g[~kept].boundary.plot(ax=ax, aspect=None, color="0.5", lw=lw)
        if style == "transcripts":
            g[kept].boundary.plot(ax=ax, aspect=None, color="white", lw=3 * lw)
        g[kept].boundary.plot(ax=ax, aspect=None, color="black", lw=1.2 * lw if style == "transcripts" else lw)
    ax.plot([xlim[1] - 30, xlim[1] - 5], [ylim[1] - 5] * 2, color="k", lw=3 * lw, solid_capstyle="butt")
    ax.set(xlim=xlim, ylim=(ylim[1], ylim[0]), aspect="equal")
    ax.axis("off")


def plot_row(sdata, sample, methods, box, style="outlines", qc="post", w=0.8, gap=0.06):
    """One row of methods at print size."""
    W, H = len(methods) * (w + gap) - gap, w + 0.3
    fig = plt.figure(figsize=(W, H))
    for i, m in enumerate(methods):
        ax = fig.add_axes([i * (w + gap) / W, 0, w / W, w / H])
        draw(ax, sdata, sample, m, box, style, qc, lw=0.2)
        ax.set_title(_constants.method_names[m].replace(" (", "\n("), fontsize=7, pad=4, linespacing=1)
    return fig


def plot_grid(sdata, sample, box, style="outlines", ncol=6):
    """All methods of a sample in one grid."""
    methods = [m for m in _constants.method_colors if f"boundaries_{m}" in sdata.shapes and _adata_path(sample, m).exists()]
    nrow = -(-len(methods) // ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(3 * ncol, 3.4 * nrow), gridspec_kw=dict(wspace=0.02, hspace=0.12))
    for ax in axs.flat:
        ax.axis("off")
    for ax, m in zip(axs.flat, methods):
        draw(ax, sdata, sample, m, box, style)
        ax.set_title(_constants.method_names[m], fontsize=14, pad=3)
    return fig
