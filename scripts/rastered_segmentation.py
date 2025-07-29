import argparse
from os.path import join

from pandas import read_csv
from sopa import aggregate, make_image_patches
from sopa.io import merscope

parser = argparse.ArgumentParser(
    description="Compute segmentation based on rasterization."
)
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument("save_path", help="Path to output folder.")
parser.add_argument("width", type=int, help="width of patches.")
parser.add_argument("overlap", type=int, help="patch overlap.")
parser.add_argument("scale", choices=["pixel", "microns"], help="Unit of measure.")
parser.add_argument("intens_rat", type=float, help="intensity ratio.")
parser.add_argument(
    "--explorer", type=bool, default=False, help="Compute explorer files."
)
args = parser.parse_args()

sdata = merscope(args.data_path)
translation = read_csv(
    join(args.data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
    sep=" ",
    header=None,
)

if args.scale == "pixel":
    make_image_patches(sdata, patch_width=args.width, patch_overlap=args.overlap)
elif args.scale == "microns":
    make_image_patches(
        sdata,
        patch_width=translation.loc[0, 0] * args.width,
        patch_overlap=translation.loc[1, 1] * args.overlap,
    )

aggregate(sdata, shapes_key="image_patches", min_intensity_ratio=args.intens_rat)
sdata["rastered_boundaries"] = sdata["image_patches"]
sdata.write(join(args.save_path, "sdata.zarr"), overwrite=True)

if args.explorer:
    from sopa.io.explorer import write

    write(
        join(args.save_path, "sdata.explorer"),
        sdata,
        shapes_key="image_patches",
        gene_column="gene",
        save_h5ad=True,
        pixel_size=1 / translation.loc[0, 0],
    )
