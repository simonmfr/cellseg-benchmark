import os
import sys
from os.path import join
from subprocess import run

import sopa
from spatialdata import read_zarr

data_path = sys.argv[1]
save_path = sys.argv[2]
staining = sys.argv[3]


def main(data_path, save_path, staining):
    """Cellpose2 algorithm by sopa with dask backend parallelized."""
    sdata = sopa.io.merscope(data_path)

    sdata.write(join(save_path, "sdata_tmp.zarr"), overwrite=True)
    sdata = read_zarr(join(save_path, "sdata_tmp.zarr"))

    sopa.make_image_patches(sdata, patch_width=8900, patch_overlap=178)

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = int(
        os.getenv("SLURM_JOB_NUM_NODES", 1)
    ) * int(os.getenv("SLURM_NTASKS_PER_NODE", 1))

    sopa.segmentation.cellpose(
        sdata,
        channels=[staining, "DAPI"],
        pretrained_model=join(
            "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/custom_models",
            f"CP_DAPI_{staining}",
        ),
        diameter=89,
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
        pixel_size=0.108,
    )

    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(join(save_path, "sdata.zarr"), overwrite=True)
    run(["rm", "-r", join(save_path, "sdata_tmp.zarr")])


if __name__ == "__main__":
    main(data_path, save_path, staining)
