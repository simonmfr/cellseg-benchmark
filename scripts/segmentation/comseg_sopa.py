#!/usr/bin/env python
import argparse
import os
import pathlib
import subprocess

import pandas as pd
import sopa
from spatialdata import read_zarr

BASE_PATH = (pathlib.Path(__file__).parents[2] / "data").resolve()

parser = argparse.ArgumentParser(description="Run ComSeg segmentation.")
parser.add_argument("data_path", help="Path to merfish output folder.")
parser.add_argument("sample", help="Sample name.")
parser.add_argument(
    "base_segmentation", help="Name of prior segmentation to use for initialization."
)
args = parser.parse_args()


def main(data_path, sample, base_segmentation):
    """ComSeg algorithm by sopa with dask backend parallelized. Resumes from cached patches when rerun."""
    path = pathlib.Path(BASE_PATH, "samples", sample, "results")
    result_dir = path / f"ComSeg_{base_segmentation}"
    tmp_path = result_dir / "sdata_tmp.zarr"

    if not tmp_path.exists():
        sdata_tmp = sopa.io.merscope(data_path)
        sdata = read_zarr(path / base_segmentation / "sdata.zarr")
        sdata[list(sdata_tmp.images.keys())[0]] = sdata_tmp[
            list(sdata_tmp.images.keys())[0]
        ]
        sdata[list(sdata_tmp.points.keys())[0]] = sdata_tmp[
            list(sdata_tmp.points.keys())[0]
        ]
        sdata.attrs["cell_segmentation_image"] = sdata_tmp.attrs[
            "cell_segmentation_image"
        ]
        sdata.attrs["transcripts_dataframe"] = sdata_tmp.attrs["transcripts_dataframe"]
        del sdata_tmp
        sdata.write(tmp_path)
    sdata = read_zarr(tmp_path)

    if "transcripts_patches" not in sdata.shapes:
        sopa.make_transcript_patches(
            sdata,
            points_key=list(sdata.points.keys())[0],
            patch_width=200,
            patch_overlap=50,
            prior_shapes_key="cellpose_boundaries",
            write_cells_centroids=True,
        )

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = int(
        os.getenv("SLURM_CPUS_PER_TASK", 1)
    )
    sopa.settings.dask_client_kwargs["threads_per_worker"] = 1

    config = str(pathlib.Path(__file__).parents[2] / "configs" / "comseg.json")
    try:
        sopa.segmentation.comseg(sdata, config=config, min_area=10, recover=True)
    except TimeoutError:
        sopa.settings.parallelization_backend = None
        sopa.segmentation.comseg(sdata, config=config, min_area=10, recover=True)

    sopa.aggregate(
        sdata,
        gene_column="gene",
        aggregate_channels=True,
        min_transcripts=10,
        points_key=list(sdata.points.keys())[0],
        image_key=list(sdata.images.keys())[0],
    )
    translation = pd.read_csv(
        pathlib.Path(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )
    sopa.io.explorer.write(
        result_dir / "sdata.explorer",
        sdata,
        points_key=list(sdata.points.keys())[0],
        image_key=list(sdata.images.keys())[0],
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )

    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(result_dir / "sdata.zarr", overwrite=True)
    subprocess.run(["rm", "-r", tmp_path])


if __name__ == "__main__":
    main(args.data_path, args.sample, args.base_segmentation)
