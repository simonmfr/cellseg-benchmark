import argparse
from os.path import join
from sopa import aggregate
from sopa.io import merscope

import sopa.aggregation.transcripts as tr
from scipy.sparse import coo_matrix as _coo

from cellseg_benchmark.sdata_utils import add_visium_boundaries

parser = argparse.ArgumentParser(
    description="Create Visium-like spot boundaries and aggregate transcripts."
)
parser.add_argument("data_path", help="Path to MERSCOPE data folder.")
parser.add_argument("save_path", help="Path to output folder.")
parser.add_argument("--points-key", default="", help="Key in sdata.points (default: first).")
parser.add_argument("--out-name", default="boundaries", help="Name for output shapes element.")
parser.add_argument("--ccd-um", type=float, default=100.0, help="Spot center-to-center distance (um).")
parser.add_argument("--diameter-um", type=float, default=55.0, help="Spot diameter (um).")
parser.add_argument("--circle-resolution", type=int, default=16, help="Circle buffer resolution.")
parser.add_argument("--explorer", action="store_true", help="Write explorer files.")
args = parser.parse_args()

sdata = merscope(args.data_path)

points_key = args.points_key or next(iter(sdata.points))
sdata = add_visium_boundaries(
    sdata,
    out_name=args.out_name,
    points_key=points_key,
    ccd_um=args.ccd_um,
    diameter_um=args.diameter_um,
    circle_resolution=args.circle_resolution,
)

# patch for older sopa versions that expect coo_matrix -> csr
if hasattr(tr, "coo_matrix"):
    def _coo_as_csr(*a, **k):
        return _coo(*a, **k).tocsr()
    tr.coo_matrix = _coo_as_csr

aggregate(sdata, shapes_key=args.out_name)
sdata.write(join(args.save_path, "sdata.zarr"), overwrite=True)

if args.explorer:
    from pandas import read_csv
    from sopa.io.explorer import write

    translation = read_csv(
        join(args.data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )

    write(
        join(args.save_path, "sdata.explorer"),
        sdata,
        gene_column="gene",
        save_h5ad=True,
        pixel_size=1 / translation.loc[0, 0],
    )
