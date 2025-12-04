import argparse
import logging
import os
import shutil
from os.path import join
from pathlib import Path
from subprocess import run

import pandas as pd
import sopa
from geopandas import GeoDataFrame
from pandas import read_csv
from shapely import Polygon
from shapely.affinity import affine_transform
from sopa._constants import SopaAttrs, SopaKeys
from sopa.aggregation.aggregation import add_standardized_table
from sopa.segmentation._transcripts import _check_transcript_patches
from sopa.segmentation.methods._proseg import (
    _read_proseg,
    _run_proseg,
    _use_zarr_output,
)
from sopa.segmentation.methods._utils import _get_executable_path
from sopa.utils import (
    delete_transcripts_patches_dirs,
    get_feature_key,
    get_transcripts_patches_dirs,
)
from spatialdata import SpatialData, read_zarr
from spatialdata.models import ShapesModel
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

log = logging.getLogger(__name__)


def _get_proseg_command(
    sdata: SpatialData, points_key: str, command_line_suffix: str
) -> str:
    proseg_executable = _get_executable_path("proseg", ".cargo")
    prior_shapes_key = sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES][
        SopaKeys.PRIOR_SHAPES_KEY
    ].iloc[0]
    use_zarr = _use_zarr_output(proseg_executable)
    feature_key = get_feature_key(sdata[points_key], raise_error=True)
    return f"proseg transcripts.csv -x x -y y -z global_z --gene-column {feature_key} --cell-id-column {prior_shapes_key} --cell-id-unassigned 0 {'--exclude-spatialdata-transcripts' if use_zarr else ''} {command_line_suffix}"


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
    assert shutil.which("proseg") is not None, (
        "Proseg is not installed. Install it according to https://github.com/dcjones/proseg"
    )

    _check_transcript_patches(sdata)

    points_key = sdata[SopaKeys.TRANSCRIPTS_PATCHES][SopaKeys.POINTS_KEY].iloc[0]

    patches_dirs = get_transcripts_patches_dirs(sdata)
    assert len(patches_dirs) == 1, (
        "Proseg is fast enough to work on a single patch. Re-run `sopa.make_transcript_patches` with `patch_width=None` and a `prior_shapes_key`."
    )
    patch_dir = Path(patches_dirs[0])

    proseg_command = _get_proseg_command(sdata, points_key, command_line_suffix)

    _run_proseg(proseg_command, patch_dir)
    adata, geo_df = _read_proseg(sdata, patch_dir, points_key)

    add_standardized_table(sdata, adata, geo_df, key_added, SopaKeys.TABLE)

    sdata.attrs[SopaAttrs.BOUNDARIES] = key_added

    if delete_cache:
        delete_transcripts_patches_dirs(sdata)

    log.info(
        "Proseg table and boundaries added (running `sopa.aggregate` is not mandatory)."
    )


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
    sdata['boundaries'] = sdata[boundaries_key]
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
        save_path / f"Proseg_3D_{base_segmentation}" / "sdata_tmp.zarr", overwrite=True
    )
    sdata = read_zarr(save_path / f"Proseg_3D_{base_segmentation}" / "sdata_tmp.zarr")

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
        sdata, patch_width=None, prior_shapes_key="auto"
    )

    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs["n_workers"] = int(
        os.getenv("SLURM_JOB_NUM_NODES", 1)
    ) * int(os.getenv("SLURM_NTASKS_PER_NODE", 1))

    proseg(sdata, delete_cache=False, command_line_suffix=proseg_flags)
    sopa.aggregate(
        sdata, gene_column="gene", aggregate_channels=True, min_transcripts=10
    )
    sopa.io.explorer.write(
        save_path / f"Proseg_3D_{base_segmentation}" / "sdata.explorer",
        sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )

    cache_dir = sopa.utils.get_cache_dir(sdata)
    del sdata[list(sdata.images.keys())[0]], sdata[list(sdata.points.keys())[0]]
    sdata.write(
        save_path / f"Proseg_3D_{base_segmentation}" / "sdata.zarr", overwrite=True
    )
    run(
        [
            "cp",
            "-r",
            cache_dir,
            str(save_path / f"Proseg_3D_{base_segmentation}" / "sdata.zarr" / str(cache_dir).split("/")[-1]),
        ]
    )
    run(["rm", "-r", join(str(save_path), f"Proseg_3D_{base_segmentation}", "sdata_tmp.zarr")])


if __name__ == "__main__":
    main(
        args.data_path,
        args.sample,
        proseg_flags,
        base_segmentation=args.base_segmentation,
    )
