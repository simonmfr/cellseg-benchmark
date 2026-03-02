import argparse
import os
from os.path import join
from pathlib import Path
from subprocess import run

import sopa
from pandas import read_csv
from spatialdata import read_zarr
from spatialdata_io import merscope

parser = argparse.ArgumentParser(
    description="Compute ProSeg segmentation without any prior segmentation."
)
parser.add_argument("data_path", help="Path to data folder.")
parser.add_argument("sample", help="Sample name.")
parser.add_argument(
    "base_segmentation", help="prior segmentation to use for initialisaton."
)
parser.add_argument(
    "proseg_flags", nargs=argparse.REMAINDER, help="Additional flags to pass to proseg."
)
args = parser.parse_args()

proseg_flags = " ".join(args.proseg_flags)


def main(data_path, sample, proseg_flags, base_segmentation):
    """Proseg 3D with vpt 3D segmentation."""
    save_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples") / sample / "results" / base_segmentation
    sdata = merscope(
        data_path,
        transcripts=True,
        mosaic_images=True,
        cells_boundaries=True,
        vpt_outputs={
            "cell_by_gene": save_path / "analysis_outputs" / "cell_by_gene.csv",
            "cell_metadata": save_path / "analysis_outputs" / "cell_metadata.csv",
            "cell_boundaries": save_path / "analysis_outputs" / "cellpose2_micron_space.parquet",
        },
    )

    image_key = list(sdata.images.keys())[0]
    boundaries_key = list(sdata.shapes.keys())[0]
    transcripts_key = list(sdata.points.keys())[0]

    sdata['transcripts'] = sdata[transcripts_key]
    sdata['cellpose_boundaries'] = sdata[boundaries_key]
    sdata['image'] = sdata[image_key]

    del sdata[image_key], sdata[boundaries_key], sdata[transcripts_key]

    sdata.attrs["cell_segmentation_image"] = 'image'
    sdata.attrs["transcripts_dataframe"] = 'transcripts'
    sdata.attrs["transcript_to_cell_assignment"] = ['cell_id', -1]
    translation = read_csv(
        join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )

    sdata.write(
        save_path / f"Proseg_3D_vpt3D_{base_segmentation}" / "sdata_tmp.zarr", overwrite=True
    )
    sdata = read_zarr(save_path / f"Proseg_3D_vpt3D_{base_segmentation}" / "sdata_tmp.zarr")

    sopa.make_transcript_patches(
        sdata, patch_width=None, prior_shapes_key="cellpose_boundaries"
    )

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = int(
        os.getenv("SLURM_JOB_NUM_NODES", 1)
    ) * int(os.getenv("SLURM_NTASKS_PER_NODE", 1))

    sopa.segmentation.proseg(
        sdata, delete_cache=False, command_line_suffix=proseg_flags
    )
    sopa.aggregate(
        sdata, gene_column="gene", aggregate_channels=True, min_transcripts=10
    )
    sopa.io.explorer.write(
        save_path / f"Proseg_3D_vpt3D_{base_segmentation}" / "sdata.explorer",
        sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )

    cache_dir = sopa.utils.get_cache_dir(sdata)
    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(
        save_path / f"Proseg_3D_vpt3D_{base_segmentation}" / "sdata.zarr", overwrite=True
    )
    run(
        [
            "cp",
            "-r",
            cache_dir,
            str(save_path / f"Proseg_3D_vpt3D_{base_segmentation}" / "sdata.zarr" / str(cache_dir).split("/")[-1]),
        ]
    )
    run(["rm", "-r", join(str(save_path), f"Proseg_3D_vpt3D_{base_segmentation}", "sdata_tmp.zarr")])


if __name__ == "__main__":
    main(
        args.data_path,
        args.sample,
        proseg_flags,
        base_segmentation=args.base_segmentation,
    )
