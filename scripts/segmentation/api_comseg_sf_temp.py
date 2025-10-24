import os
import sys
import time
from os.path import join
from subprocess import run

import sopa
from pandas import read_csv
from spatialdata import read_zarr

data_path = sys.argv[1]
sample = sys.argv[2]
base_segmentation = sys.argv[3]


def main(data_path, sample, base_segmentation):
    """ComSeg algorithm by sopa with dask backend parallelized."""
    print("Reading raw MERSCOPE data...")
    start = time.time()
    sdata_tmp = sopa.io.merscope(data_path)  # to read in the images and points
    print(f"Reading done in {time.time() - start:.1f}s")

    path = f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{sample}/results"

    print("Loading base segmentation...")
    sdata = read_zarr(
        join(path, base_segmentation, "sdata.zarr")
    )  # enthält keine Bilder oder transcripte
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
    # backing für memory efficiency
    start = time.time()
    sdata.write(
        join(path, f"ComSeg_{base_segmentation}", "sdata_tmp.zarr"), overwrite=True
    )
    print(f"Writing done in {time.time() - start:.1f}s")
    sdata = read_zarr(join(path, f"ComSeg_{base_segmentation}", "sdata_tmp.zarr"))

    print("Creating transcript patches using Cellpose boundaries...")
    start = time.time()
    # Annahme: nur cellpose prior wird benutzt
    sopa.make_transcript_patches(
        sdata,
        patch_width=200,
        patch_overlap=50,
        prior_shapes_key="cellpose_boundaries",
        write_cells_centroids=True,
    )
    print(f"Creating transcript patches done in {time.time() - start:.1f}s")

    sopa.settings.parallelization_backend = "dask"
    n_cpus = int(os.getenv("SLURM_CPUS_PER_TASK", 1))
    sopa.settings.dask_client_kwargs["n_workers"] = n_cpus
    sopa.settings.dask_client_kwargs["threads_per_worker"] = 1
    # sopa.settings.dask_client_kwargs["n_workers"] = int(
    #    os.getenv("SLURM_JOB_NUM_NODES", 1)
    # ) * int(os.getenv("SLURM_NTASKS_PER_NODE", 1))
    sopa.settings.dask_client_kwargs["timeout"] = "600000"

    print("Running ComSeg segmentation...")
    path_json = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/comseg.json"
    start = time.time()
    sopa.segmentation.comseg(sdata, config=path_json, min_area=10, delete_cache=False)
    print(f"Segmentation done in {time.time() - start:.1f}s")

    print("Aggregating transcript data...")
    start = time.time()
    sopa.aggregate(
        sdata, gene_column="gene", aggregate_channels=True, min_transcripts=10
    )
    print(f"Aggregation done in {time.time() - start:.1f}s")

    print("Writing Sopa Explorer output...")
    start = time.time()
    sopa.io.explorer.write(
        join(path, f"ComSeg_{base_segmentation}", "sdata.explorer"),
        sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )
    print(f"Explorer write done in {time.time() - start:.1f}s")

    print("Cleaning up temporary data and saving final result...")
    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(join(path, f"ComSeg_{base_segmentation}", "sdata.zarr"))
    run(["rm", "-r", join(path, f"ComSeg_{base_segmentation}", "sdata_tmp.zarr")])

    print("Done.")


if __name__ == "__main__":
    main(data_path, sample, base_segmentation)
