import argparse
import os
from os.path import join
from subprocess import run

import sopa
from pandas import read_csv
from spatialdata import read_zarr

parser = argparse.ArgumentParser(description="Compute ComSeg segmentation.")
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument("sample", help="sample name.")
parser.add_argument(
    "base_segmentation", help="prior segmentation to use for initialisaton."
)
args = parser.parse_args()


def main(data_path, sample, base_segmentation):
    """ComSeg algorithm by sopa with dask backend parallelized."""
    sdata_tmp = sopa.io.merscope(data_path)  # to read in the images and points
    path = f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{sample}/results"
    sdata = read_zarr(
        join(path, base_segmentation, "sdata.zarr")
    )  # enthält keine Bilder oder transcripte
    sdata[list(sdata_tmp.images.keys())[0]] = sdata_tmp[
        list(sdata_tmp.images.keys())[0]
    ]
    sdata[list(sdata_tmp.points.keys())[0]] = sdata_tmp[
        list(sdata_tmp.points.keys())[0]
    ]
    translation = read_csv(
        join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )
    del sdata_tmp

    # backing für memory efficiency
    sdata.write(
        join(path, f"ComSeg_{base_segmentation}", "sdata_tmp.zarr"), overwrite=True
    )
    sdata = read_zarr(join(path, f"ComSeg_{base_segmentation}", "sdata_tmp.zarr"))

    # Annahme: nur cellpose prior wird benutzt
    sopa.make_transcript_patches(
        sdata,
        patch_width=1200,
        patch_overlap=50,
        prior_shapes_key="cellpose_boundaries",
        write_cells_centroids=True,
    )

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = int(
        os.getenv("SLURM_JOB_NUM_NODES", 1)
    ) * int(os.getenv("SLURM_NTASKS_PER_NODE", 1))
    sopa.settings.dask_client_kwargs["timeout"] = "600000"

    path_json = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/comseg.json"
    sopa.segmentation.comseg(sdata, config=path_json, min_area=10, delete_cache=False)

    sopa.aggregate(
        sdata, gene_column="gene", aggregate_channels=True, min_transcripts=10
    )
    sopa.io.explorer.write(
        join(path, f"ComSeg_{base_segmentation}", "sdata.explorer"),
        sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )

    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(join(path, f"ComSeg_{base_segmentation}", "sdata.zarr"))
    run(["rm", "-r", join(path, f"ComSeg_{base_segmentation}", "sdata_tmp.zarr")])


if __name__ == "__main__":
    main(args.data_path, args.sample, args.base_segmentation)
