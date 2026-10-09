#!/usr/bin/env python
import argparse

import spatialdata as sd

from cellseg_benchmark.sdata_utils import assign_transformations

parser = argparse.ArgumentParser(
    description="Rewrite transformations of a master sdata: 2D images/shapes get 2D transformations, global = pixel."
)
parser.add_argument("sdata_path", help="Path to master sdata, e.g. samples/<sample>/sdata_z3.zarr.")
args = parser.parse_args()

sdata = sd.read_zarr(args.sdata_path, selection=("images", "points", "shapes"))
micron_to_pixel = sd.transformations.Affine(
    sd.transformations.get_transformation(sdata[list(sdata.points.keys())[0]], "pixel").to_affine_matrix(("x", "y"), ("x", "y")),
    ("x", "y"),
    ("x", "y"),
)

for key, image in sdata.images.items():
    sd.transformations.set_transformation(
        image,
        {
            "micron": micron_to_pixel.inverse(),
            "pixel": sd.transformations.Identity(),
            "global": sd.transformations.Identity(),
        },
        set_all=True,
    )
    sdata.write_transformations(key)

for key, boundaries in sdata.shapes.items():
    assign_transformations(boundaries, key.removeprefix("boundaries_"), micron_to_pixel)
    sdata.write_transformations(key)
    print(f"{key}: done")
