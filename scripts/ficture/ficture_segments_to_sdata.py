#!/usr/bin/env python
"""Turn a FICTURE pixel decode into cell boundaries as SpatialData objects.

For one sample two boundary sets are written, each as ``shapes["boundaries"]``
in a ``sdata.zarr`` (polygons in micrometers, Identity transform to "micron"):

* ``Ficture_segments``       - raw FICTURE segments (one polygon per connected,
  same-factor patch; neighbouring cells of one type stay fused).
* ``Ficture_segments_dapi``  - the same segments split into single cells using
  DAPI nuclei as seeds. A FICTURE segment carries no membrane, so its shape is
  not a reliable border; instead every pixel is assigned to its nearest nucleus
  within the segment. A segment with N nuclei becomes N cells; a segment with
  one nucleus stays whole; nucleus-free segments are kept (n_nuclei == 0).

Logic lives in cellseg_benchmark.ficture_utils; this is the per-sample CLI:

    python ficture_segments_to_sdata.py <SAMPLE_ID>
"""
import argparse
import logging
import os
from pathlib import Path

import dask

dask.config.set({"dataframe.query-planning": False})  # required before spatialdata import

import anndata as ad  # noqa: E402

from cellseg_benchmark import _constants  # noqa: E402
from cellseg_benchmark.ficture_utils import (  # noqa: E402
    aggregate_tables,
    build_factor_raster,
    plot_qc,
    segments_to_boundaries,
    split_by_nuclei,
    write_boundaries,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("ficture_segments")


def main():
    """Build FICTURE boundary sdatas (and transcript tables) for one sample."""
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("sample", help="sample id, e.g. foxf2_s7_r1")
    p.add_argument("--base", default=_constants.BASE_PATH, help="cellseg-benchmark base path")
    p.add_argument("--data-path", default=None,
                   help="MERSCOPE output folder; if given, aggregate transcripts into a table")
    p.add_argument("--min-transcripts", type=int, default=10,
                   help="drop cells with fewer transcripts when building the table")
    p.add_argument("--res", type=float, default=1.5, help="grid size in um")
    p.add_argument("--min-um2", type=float, default=5.0, help="drop blobs / fill holes below this")
    p.add_argument("--connectivity", type=int, default=1, choices=(1, 2),
                   help="1 = 4-conn (avoids diagonal 1px bridges), 2 = 8-conn")
    p.add_argument("--n-jobs", type=int,
                   default=int(os.environ.get("SLURM_CPUS_PER_TASK")
                               or os.environ.get("SLURM_CPUS_ON_NODE")
                               or os.cpu_count() or 8),
                   help="parallel workers (default: allocated SLURM cpus, else all cores)")
    p.add_argument("--no-plot", action="store_true", help="skip QC plots")
    args = p.parse_args()

    root = Path(args.base) / "samples" / args.sample / "results"
    pix = root / "Ficture" / "output" / "decode.pixel.sorted.tsv.gz"
    dapi = root / "vpt_3D_DAPI_nuclei" / "sdata.zarr" / "tables" / "table"
    out_raw = root / "Ficture_segments"
    out_split = root / "Ficture_segments_dapi"

    log.info("[%s] rasterizing FICTURE pixels (res=%s um)", args.sample, args.res)
    lab, affine, bbox = build_factor_raster(str(pix), args.res, args.min_um2)
    aspect = (bbox[3] - bbox[2]) / (bbox[1] - bbox[0])

    # raw FICTURE segments
    seg = segments_to_boundaries(lab, affine)
    log.info("[%s] %d unsplit segments", args.sample, len(seg))

    # nuclei-split cells
    log.info("[%s] loading nuclei and splitting", args.sample)
    obs = ad.read_zarr(str(dapi)).obs
    nuclei = obs[["center_x", "center_y"]].to_numpy(float)
    cells = split_by_nuclei(lab, affine, nuclei, obs.index.to_numpy(),
                            args.connectivity, args.n_jobs)
    log.info("[%s] %d cells (%d nucleated, %d nucleus-free)", args.sample, len(cells),
             int((cells.n_nuclei == 1).sum()), int((cells.n_nuclei == 0).sum()))

    # write boundaries (+ transcript table if --data-path given)
    out_raw.mkdir(parents=True, exist_ok=True)
    out_split.mkdir(parents=True, exist_ok=True)
    if args.data_path:
        log.info("[%s] aggregating transcripts into tables", args.sample)
        aggregate_tables(args.data_path,
                         [(seg, out_raw / "sdata.zarr"), (cells, out_split / "sdata.zarr")],
                         min_transcripts=args.min_transcripts, n_workers=args.n_jobs)
    else:
        write_boundaries(seg, out_raw / "sdata.zarr")
        write_boundaries(cells, out_split / "sdata.zarr")

    # QC plots
    if not args.no_plot:
        plot_qc(seg, None, aspect, f"{args.sample}: FICTURE segments (unsplit)",
                out_raw / "ficture_segments_qc.png")
        plot_qc(cells, nuclei, aspect, f"{args.sample}: nuclei-split cells",
                out_split / "ficture_segments_qc.png")
    log.info("[%s] done", args.sample)


if __name__ == "__main__":
    main()
