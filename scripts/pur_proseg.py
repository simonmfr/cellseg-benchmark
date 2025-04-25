import os
import sys
from os.path import join
from subprocess import run
import logging

import sopa
from pandas import read_csv
from spatialdata import read_zarr, SpatialData

import shutil
from pathlib import Path

from sopa._constants import SopaAttrs, SopaKeys
from sopa.aggregation.aggregation import add_standardized_table
from sopa.utils import (
    delete_transcripts_patches_dirs,
    get_feature_key,
    get_transcripts_patches_dirs,
)
from sopa.segmentation._transcripts import _check_transcript_patches
from sopa.segmentation.methods._proseg import _run_proseg, _read_proseg

data_path = sys.argv[1]
sample = sys.argv[2]
proseg_flags = " ".join(sys.argv[3:])

log = logging.getLogger(__name__)

def _get_proseg_command(sdata: SpatialData, points_key: str, command_line_suffix: str) -> str:
    feature_key = get_feature_key(sdata[points_key], raise_error=True)
    return f"proseg transcripts.csv -x x -y y -z z --gene-column {feature_key} --cell-id-column cell_id --cell-id-unassigned '-- -1' {command_line_suffix}"

def proseg(
    sdata: SpatialData,
    delete_cache: bool = True,
    command_line_suffix: str = "",
    key_added: str = SopaKeys.PROSEG_BOUNDARIES,
):
    """Run [`proseg`](https://github.com/dcjones/proseg) segmentation on a SpatialData object, and add the corresponding cell boundaries and `AnnData` table with counts.

    !!! warning "Proseg installation"
        Make sure to install [`proseg`](https://github.com/dcjones/proseg) separately before running this function.

    !!! info "Proseg usage specificities"
        Contrary to most other segmentation tools, `proseg` will only run on one patch. I.e., you need
        to run [`sopa.make_transcript_patches`](../patches/#sopa.make_transcript_patches) with `patch_width=None` and a `prior_shapes_key` before running `proseg`.

        Also, note that aggregation is not necessary after running `proseg`.

    Args:
        sdata: A `SpatialData` object.
        delete_cache: Whether to delete the cache after segmentation.
        command_line_suffix: Optional suffix to add to the proseg command line.
        key_added: Name of the shapes element to be added to `sdata.shapes`.
    """
    assert (
        shutil.which("proseg") is not None
    ), "Proseg is not installed. Install it according to https://github.com/dcjones/proseg"

    _check_transcript_patches(sdata)

    points_key = sdata[SopaKeys.TRANSCRIPTS_PATCHES][SopaKeys.POINTS_KEY].iloc[0]

    patches_dirs = get_transcripts_patches_dirs(sdata)
    assert (
        len(patches_dirs) == 1
    ), "Proseg is fast enough to work on a single patch. Re-run `sopa.make_transcript_patches` with `patch_width=None` and a `prior_shapes_key`."
    patch_dir = Path(patches_dirs[0])

    proseg_command = _get_proseg_command(sdata, points_key, command_line_suffix)

    _run_proseg(proseg_command, patch_dir)
    adata, geo_df = _read_proseg(sdata, patch_dir, points_key)

    add_standardized_table(sdata, adata, geo_df, key_added, SopaKeys.TABLE)

    sdata.attrs[SopaAttrs.BOUNDARIES] = key_added

    if delete_cache:
        delete_transcripts_patches_dirs(sdata)

    log.info("Proseg table and boundaries added (running `sopa.aggregate` is not mandatory).")

def main(data_path, sample, proseg_flags):
    """ComSeg algorithm by sopa with dask backend parallelized."""
    sdata = sopa.io.merscope(data_path)  # to read in the images and points
    path = f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{sample}/results"
    translation = read_csv(
        join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )

    # backing für memory efficiency
    sdata.write(
        join(path, f"Proseg_pure", "sdata_tmp.zarr"), overwrite=True
    )
    sdata = read_zarr(join(path, f"Proseg_pure", "sdata_tmp.zarr"))

    sopa.make_transcript_patches(sdata, patch_width=None)

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = int(
        os.getenv("SLURM_JOB_NUM_NODES", 1)
    ) * int(os.getenv("SLURM_NTASKS_PER_NODE", 1))

    proseg(
        sdata, delete_cache=False, command_line_suffix=proseg_flags
    )
    sopa.aggregate(
        sdata, gene_column="gene", aggregate_channels=True, min_transcripts=10
    )
    sopa.io.explorer.write(
        join(path, f"Proseg_pure", "sdata.explorer"),
        sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )

    cache_dir = sopa.utils.get_cache_dir(sdata)
    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(join(path, f"Proseg_pure", "sdata.zarr"))
    run(
        [
            "cp",
            "-r",
            cache_dir,
            join(
                path,
                f"Proseg_pure",
                "sdata.zarr",
                str(cache_dir).split("/")[-1],
            ),
        ]
    )
    run(["rm", "-r", join(path, f"Proseg_pure", "sdata_tmp.zarr")])


if __name__ == "__main__":
    main(data_path, sample, proseg_flags)