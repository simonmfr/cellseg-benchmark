import argparse
import os
from os.path import join
from subprocess import run

import pandas as pd
import sopa
import toml
from geopandas import GeoDataFrame
from pandas import read_csv
from shapely import Polygon
from shapely.affinity import affine_transform
from spatialdata import read_zarr
from spatialdata.models import ShapesModel

parser = argparse.ArgumentParser(description="Compute Baysor segmentation.")
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument(
    "base_segmentation", help="prior segmentation to use for initialisaton."
)
parser.add_argument("confidence", type=float, help="confidence of prior segmentation.")
parser.add_argument("sample", help="sample name.")
args = parser.parse_args()


def main(data_path, base_segmentation, confidence, sample):
    """Baysor algorithm by sopa with dask backend parallelized."""
    sdata_tmp = sopa.io.merscope(data_path)
    path = f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{sample}/results"
    sdata = read_zarr(join(path, base_segmentation, "sdata.zarr"))
    translation = read_csv(
        join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )
    sdata[list(sdata_tmp.images.keys())[0]] = sdata_tmp[
        list(sdata_tmp.images.keys())[0]
    ]
    sdata[list(sdata_tmp.points.keys())[0]] = sdata_tmp[
        list(sdata_tmp.points.keys())[0]
    ]
    del sdata_tmp

    # backing für memory efficiency
    sdata.write(
        join(path, f"Baysor_2D_{base_segmentation}_{confidence}", "sdata_tmp.zarr"),
        overwrite=True,
    )
    sdata = read_zarr(
        join(path, f"Baysor_2D_{base_segmentation}_{confidence}", "sdata_tmp.zarr")
    )

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

    sopa.make_transcript_patches(
        sdata,
        patch_width=1000,
        patch_overlap=20,
        prior_shapes_key="cellpose_boundaries",
        points_key=list(sdata.points.keys())[0],
    )

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = int(
        os.getenv("SLURM_JOB_NUM_NODES", 1)
    ) * int(os.getenv("SLURM_NTASKS_PER_NODE", 1))

    path_toml = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/baysor_2D_config.toml"
    with open(path_toml, "r") as f:
        config = toml.load(f)
    config["segmentation"]["prior_segmentation_confidence"] = confidence
    sopa.segmentation.baysor(sdata, config=config, delete_cache=True, force=True)

    sopa.aggregate(
        sdata,
        gene_column="gene",
        aggregate_channels=True,
        min_transcripts=10,
        points_key=list(sdata.points.keys())[0],
        image_key=list(sdata.images.keys())[0],
    )
    sopa.io.explorer.write(
        join(path, f"Baysor_2D_{base_segmentation}_{confidence}", "sdata.explorer"),
        sdata,
        points_key=list(sdata.points.keys())[0],
        image_key=list(sdata.images.keys())[0],
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=0.108,
    )

    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(
        join(path, f"Baysor_2D_{base_segmentation}_{confidence}", "sdata.zarr"),
        overwrite=True,
    )
    run(
        [
            "rm",
            "-r",
            join(path, f"Baysor_2D_{base_segmentation}_{confidence}", "sdata_tmp.zarr"),
        ]
    )


if __name__ == "__main__":
    main(args.data_path, args.base_segmentation, args.confidence, args.sample)
