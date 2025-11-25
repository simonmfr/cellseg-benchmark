import argparse
import time
from os.path import join
from subprocess import run

import sopa
from pandas import read_csv
from spatialdata import read_zarr

parser = argparse.ArgumentParser(description="Run ComSeg segmentation.")
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument("sample", help="Sample name.")
parser.add_argument(
    "base_segmentation", help="Name of prior segmentation to use for initializaton."
)
args = parser.parse_args()


def main(data_path, sample, base_segmentation):
    """ComSeg algorithm by sopa with dask backend parallelized."""
    print("Reading raw MERSCOPE data...")
    sdata_tmp = sopa.io.merscope(data_path)  # to read in the images and points
    path = f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{sample}/results"

    print("Loading base segmentation...")
    sdata = read_zarr(join(path, base_segmentation, "sdata.zarr"))
    sdata[list(sdata_tmp.images.keys())[0]] = sdata_tmp[
        list(sdata_tmp.images.keys())[0]
    ]
    sdata[list(sdata_tmp.points.keys())[0]] = sdata_tmp[
        list(sdata_tmp.points.keys())[0]
    ]

    # fix legacy path inconsistency
    sdata.attrs["cell_segmentation_image"] = (
        "_".join(data_path.rstrip("/").split("/")[-2:]) + "_z3"
    )
    sdata.attrs["transcripts_dataframe"] = (
        "_".join(data_path.rstrip("/").split("/")[-2:]) + "_transcripts"
    )

    translation = read_csv(
        join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )
    del sdata_tmp

    print("Backing sdata to temporary Zarr store...")
    start = time.time()
    sdata.write(
        join(path, f"ComSeg_{base_segmentation}", "sdata_tmp.zarr"), overwrite=True
    )
    print(f"Writing done in {(time.time() - start) / 3600:.2f}h")
    sdata = read_zarr(join(path, f"ComSeg_{base_segmentation}", "sdata_tmp.zarr"))

    print("Creating transcript patches using Cellpose boundaries...")
    start = time.time()
    sopa.make_transcript_patches(
        sdata,
        patch_width=1000,
        patch_overlap=20,
        prior_shapes_key="cellpose_boundaries",
        write_cells_centroids=True,
    )
    print(f"Creating transcript patches done in {(time.time() - start) / 3600:.2f}h")

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = 5
    sopa.settings.dask_client_kwargs["threads_per_worker"] = 4
    sopa.settings.dask_client_kwargs["timeout"] = "600000"

    # sopa.settings.parallelization_backend = "dask"
    # sopa.settings.dask_client_kwargs["n_workers"] = int(
    #    os.getenv("SLURM_JOB_NUM_NODES", 1)
    # ) * int(os.getenv("SLURM_NTASKS_PER_NODE", 1))
    # sopa.settings.dask_client_kwargs["timeout"] = "600000"

    print("Running ComSeg segmentation...")
    path_json = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/comseg.json"
    start = time.time()
    try:
        sopa.segmentation.comseg(
            sdata, config=path_json, min_area=10, delete_cache=False
        )
    except Exception as e:
        print("WARNING: comseg finished but Dask cleanup threw:", repr(e))
    print(f"Segmentation done in {(time.time() - start) / 3600:.2f}h")

    print("Aggregating transcript data...")
    start = time.time()
    sopa.aggregate(
        sdata, gene_column="gene", aggregate_channels=True, min_transcripts=10
    )
    print(f"Aggregation done in {(time.time() - start) / 3600:.2f}h")

    print("Writing Sopa Explorer output...")
    start = time.time()
    sopa.io.explorer.write(
        join(path, f"ComSeg_{base_segmentation}", "sdata.explorer"),
        sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )
    print(f"Explorer write done in {(time.time() - start) / 3600:.2f}h")

    print("Cleaning up temporary data and saving final result...")
    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(join(path, f"ComSeg_{base_segmentation}", "sdata.zarr"))
    run(["rm", "-r", join(path, f"ComSeg_{base_segmentation}", "sdata_tmp.zarr")])

    print("Done.")


if __name__ == "__main__":
    main(args.data_path, args.sample, args.base_segmentation)
