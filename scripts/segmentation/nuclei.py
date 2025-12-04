import argparse
import os
from os.path import join
from subprocess import run

import pandas as pd
import sopa
from geopandas import GeoDataFrame
from pandas import read_csv
from shapely import Polygon
from shapely.affinity import affine_transform
from spatialdata import read_zarr
from spatialdata.models import ShapesModel

parser = argparse.ArgumentParser(description="Compute Cellpose nuclei segmentation.")
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument("save_path", help="Path to output folder.")
args = parser.parse_args()


def main(data_path, save_path):
    """Cellpose nuclei algorithm by sopa with dask backend parallelized."""
    sdata = sopa.io.merscope(data_path)
    translation = read_csv(
        join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )

    sdata.write(join(save_path, "sdata_tmp.zarr"), overwrite=True)
    sdata = read_zarr(join(save_path, "sdata_tmp.zarr"))

    if "ABCAtlas" in data_path:
        coords = pd.read_csv(
            join("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/ABC_explorers/",
                 f"{data_path.split('/')[-1]}_ROI.csv"),
            skiprows=2
            )
        polygon = Polygon([
            (x, y) for x, y in coords.values
        ])
        polygon_spat = affine_transform(polygon,
                                        [translation.iloc[0, 0], translation.iloc[0, 1], translation.iloc[1, 0],
                                         translation.iloc[1, 1], translation.iloc[0, 2], translation.iloc[1, 2]])
        gdf = GeoDataFrame({'geometry': [polygon_spat]}, geometry='geometry')
        sdata['region_of_interest'] = ShapesModel(gdf)

    sopa.make_image_patches(sdata, patch_width=8900, patch_overlap=178)

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = int(
        os.getenv("SLURM_JOB_NUM_NODES", 1)
    ) * int(os.getenv("SLURM_NTASKS_PER_NODE", 1))

    sopa.segmentation.cellpose(
        sdata,
        channels=["DAPI"],
        model_type="nuclei",
        diameter=60,
        flow_threshold=2,
        cellprob_threshold=-6,
        min_area=2000,
    )

    sopa.aggregate(
        sdata, gene_column="gene", aggregate_channels=True, min_transcripts=10
    )
    sopa.io.explorer.write(
        join(save_path, "sdata.explorer"),
        sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )

    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(join(save_path, "sdata.zarr"), overwrite=True)
    run(["rm", "-r", join(save_path, "sdata_tmp.zarr")])


if __name__ == "__main__":
    main(args.data_path, args.save_path)
