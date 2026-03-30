import argparse
from os.path import join
from pathlib import Path
from subprocess import run

import sopa
from pandas import read_csv
from spatialdata import read_zarr

parser = argparse.ArgumentParser(description="Remove overlaps from Baysor segmentation .")
parser.add_argument("data_path", help="Path to merfish output folder.")
parser.add_argument("zarr_path", help="Path to Baysor .zarr-folder.")
parser.add_argument("save_path", help="Path to output folder.")
args = parser.parse_args()

def main(data_path, zarr_path, save_path):
    sdata_tmp = sopa.io.merscope(Path(data_path))
    sdata = read_zarr(Path(zarr_path) / "sdata.zarr")

    sdata[list(sdata_tmp.images.keys())[0]] = sdata_tmp[
        list(sdata_tmp.images.keys())[0]
    ]
    sdata[list(sdata_tmp.points.keys())[0]] = sdata_tmp[
        list(sdata_tmp.points.keys())[0]
    ]
    sdata.attrs["cell_segmentation_image"] = sdata_tmp.attrs["cell_segmentation_image"]
    sdata.attrs["transcripts_dataframe"] = sdata_tmp.attrs["transcripts_dataframe"]
    del sdata_tmp

    del sdata['table'] #otherwise transcript counts are not counted again
    boundaries = sdata['baysor_boundaries'].copy()
    boundaries = sopa.shapes.remove_overlap(boundaries) #Very slow

    sdata['baysor_boundaries'] = boundaries
    sopa.aggregate(
        sdata,
        gene_column="gene",
        aggregate_channels=True,
        min_transcripts=10,
        points_key=list(sdata.points.keys())[0],
        image_key=list(sdata.images.keys())[0],
    )

    sdata.write(Path(save_path) / "sdata.zarr", overwrite=True)
    sdata = read_zarr(Path(save_path) / "sdata.zarr")

    translation = read_csv(
            join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
            sep=" ",
            header=None,
        )
    sopa.explorer.write(
        Path(save_path) / "sdata.explorer",
        sdata,
        gene_column="gene",
        ram_threshold_gb=4,
        pixel_size=1 / translation.iloc[0, 0],
    )
    run(["rm", "-r", join(save_path, "sdata.zarr", "images")])
    run(["rm", "-r", join(save_path, "sdata.zarr", "points")])

if __name__ == "__main__":
    main(args.data_path, args.zarr_path, args.save_path)