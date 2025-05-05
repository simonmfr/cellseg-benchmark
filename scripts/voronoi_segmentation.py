import sys
from os.path import join
from pathlib import Path

import shapely
from geopandas import GeoDataFrame
from numpy.random import random
from pandas import read_csv
from scipy.spatial import Voronoi
from sopa import aggregate
from sopa.io import merscope
from sopa.io.explorer import write
from spatialdata import read_zarr
from spatialdata.models import ShapesModel
from tifffile import imread

data_path = sys.argv[1]
save_path = sys.argv[2]

sdata = merscope(data_path)
shape = imread(join(data_path, "images/mosaic_DAPI_z3.tif")).shape
translation = read_csv(
    join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
    sep=" ",
    header=None,
)

N = read_zarr(
    join(
        str(Path(save_path).parent.resolve()),
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
sdata.write(join(save_path, "sdata.zarr"), overwrite=True)

#write(
#    join(save_path, "sdata.explorer"),
#    sdata,
#    gene_column="gene",
#    save_h5ad=True,
#    pixel_size=1 / translation.loc[0, 0],
#)
