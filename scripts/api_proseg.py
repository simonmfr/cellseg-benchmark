import os
import sys
from os.path import join
from subprocess import run

import sopa
from pandas import read_csv
from spatialdata import read_zarr

data_path = sys.argv[1]
sample = sys.argv[2]
base_segmentation = sys.argv[3]
proseg_flags = " ".join(sys.argv[4:])


def main(data_path, sample, base_segmentation, proseg_flags):
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
        join(path, f"Proseg_{base_segmentation}", "sdata_tmp.zarr"), overwrite=True
    )
    sdata = read_zarr(join(path, f"Proseg_{base_segmentation}", "sdata_tmp.zarr"))

    # Annahme: nur cellpose prior wird benutzt
    sopa.make_transcript_patches(
        sdata,
        points_key=list(sdata.points.keys())[0],
        patch_width=None,
        prior_shapes_key="cellpose_boundaries",
        write_cells_centroids=True,
    )

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = int(
        os.getenv("SLURM_JOB_NUM_NODES", 1)
    ) * int(os.getenv("SLURM_NTASKS_PER_NODE", 1))

    sopa.segmentation.proseg(
        sdata, delete_cache=False, command_line_suffix=proseg_flags
    )
    sopa.aggregate(
        sdata, gene_column="gene", aggregate_channels=True, min_transcripts=10,
        points_key=list(sdata.points.keys())[0], image_key=list(sdata.images.keys())[0]
    )
    sopa.io.explorer.write(
        join(path, f"Proseg_{base_segmentation}", "sdata.explorer"),
        sdata,
        points_key=list(sdata.points.keys())[0],
        image_key=list(sdata.images.keys())[0],
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )

    cache_dir = sopa.utils.get_cache_dir(sdata)
    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(join(path, f"Proseg_{base_segmentation}", "sdata.zarr"), overwrite=True)
    run(
        [
            "cp",
            "-r",
            cache_dir,
            join(
                path,
                f"Proseg_{base_segmentation}",
                "sdata.zarr",
                str(cache_dir).split("/")[-1],
            ),
        ]
    )
    run(["rm", "-r", join(path, f"Proseg_{base_segmentation}", "sdata_tmp.zarr")])


if __name__ == "__main__":
    main(data_path, sample, base_segmentation, proseg_flags)
