import sys
from os.path import join

from pandas import read_csv
from sopa import aggregate, make_image_patches
from sopa.io import merscope

data_path = sys.argv[1]
save_path = sys.argv[2]
width = int(sys.argv[3])
overlap = int(sys.argv[4])
intens_rat = float(sys.argv[6])

sdata = merscope(data_path)
translation = read_csv(
    join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
    sep=" ",
    header=None,
)

if sys.argv[5] == "pixel":
    make_image_patches(sdata, patch_width=width, patch_overlap=overlap)
elif sys.argv[5] == "microns":
    make_image_patches(
        sdata,
        patch_width=translation.loc[0, 0] * width,
        patch_overlap=translation.loc[1, 1] * overlap,
    )

aggregate(sdata, shapes_key="image_patches", min_intensity_ratio=intens_rat)
sdata["rastered_boundaries"] = sdata["image_patches"]
sdata.write(join(save_path, "sdata.zarr"), overwrite=True)

# write(
#    join(save_path, "sdata.explorer"),
#    sdata,
#    shapes_key="image_patches",
#    gene_column="gene",
#    save_h5ad=True,
#    pixel_size=1 / translation.loc[0, 0],
# )
