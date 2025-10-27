import argparse
from os.path import join
from pathlib import Path

import shapely
from geopandas import GeoDataFrame
from numpy.random import random
from pandas import read_csv
from scipy.spatial import Voronoi
from sopa import aggregate
from sopa.io import merscope
from spatialdata import read_zarr
from spatialdata.models import ShapesModel
from tifffile import imread

parser = argparse.ArgumentParser(
    description="Compute Voronoi segmentation based on Negative Control 10um."
)
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument("save_path", help="Path to output folder.")
parser.add_argument(
    "--explorer", type=bool, default=False, help="Compute explorer files."
)
args = parser.parse_args()

sdata = merscope(args.data_path)
shape = imread(join(args.data_path, "images/mosaic_DAPI_z3.tif")).shape
translation = read_csv(
    join(args.data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
    sep=" ",
    header=None,
)

N = read_zarr(
    join(
        str(Path(args.save_path).parent.resolve()),
        "Negative_Control_Rastered_10",
        "sdata.zarr",
    )
)["table"].n_obs
points = random((N, 2))
points[:, 1] = (points[:, 1] * shape[0]).astype(int)
points[:, 0] = (points[:, 0] * shape[1]).astype(int)

vor = Voronoi(points)
lines = [
    shapely.geometry.LineString(vor.vertices[line])
    for line in vor.ridge_vertices
    if -1 not in line
]
polygons = shapely.ops.polygonize(lines)
gdf = GeoDataFrame({"geometry": [polygon for polygon in polygons]})

sdata["cellpose_boundaries"] = ShapesModel.parse(gdf)
aggregate(sdata, shapes_key="cellpose_boundaries")
sdata.write(join(args.save_path, "sdata.zarr"), overwrite=True)

if args.explorer:
    from sopa.io.explorer import write

    write(
        join(args.save_path, "sdata.explorer"),
        sdata,
        gene_column="gene",
        save_h5ad=True,
        pixel_size=1 / translation.loc[0, 0],
    )
