import argparse
import os
import pathlib
import subprocess
import pandas as pd
import spatialdata
import spatialdata_io
import sopa.aggregation
import sopa.io.explorer
import sopa.utils

def main():
    parser = argparse.ArgumentParser(
        description="Convert vpt 2D output to sdata."
    )
    parser.add_argument("data_path", help="Path to vpt folder.")
    parser.add_argument("save_path", help="Path to output folder.")
    parser.add_argument(
        "--explorer", action="store_true", help="if to compute explorer files"
    )
    args = parser.parse_args()

    save_path = pathlib.Path(args.save_path)
    data_path = pathlib.Path(args.data_path)

    assert any(
        "cell_by_gene.csv" in file
        for file in os.listdir(save_path / "analysis_outputs")
    ), "not correctly computed"

    sdata = spatialdata_io.merscope(
        data_path,
        vpt_outputs={
            "cell_by_gene": save_path / "analysis_outputs" / "cell_by_gene.csv",
            "cell_metadata": save_path / "analysis_outputs" / "cell_metadata.csv",
            "cell_boundaries": save_path / "analysis_outputs" / "cellpose2_micron_space.parquet",
        },
    )

    sdata["table"].obs.rename(columns={"EntityID": "cell_id"}, inplace=True)

    shapes_key = list(sdata.shapes.keys())[0]

    translation = pd.read_csv(
        data_path / "images" / "micron_to_mosaic_pixel_transform.csv",
        sep=" ",
        header=None,
    )

    boundaries = sdata[shapes_key]
    boundaries.rename(columns={"EntityID": "cell_id"}, inplace=True)
    boundaries.set_index("cell_id", drop=False, inplace=True)
    boundaries.index = boundaries.index.rename(None)

    sdata["table"].uns["spatialdata_attrs"]["instance_key"] = "cell_id"

    sdata["table"].obsm['intensities'] = pd.DataFrame(
        sopa.aggregation.aggregate_channels(sdata, shapes_key=shapes_key),
        columns=sopa.utils.validated_channel_names(
            sopa.utils.get_spatial_image(sdata, list(sdata.images.keys())[0], return_key=True)[1]
        ),
        index=sdata[shapes_key].index.astype(str)
    )

    zarr_path = save_path / "sdata.zarr"
    sdata.write(str(zarr_path), overwrite=True)
    sdata = spatialdata.read_zarr(str(zarr_path))

    if args.explorer:
        sopa.io.explorer.write(
            str(save_path / "sdata.explorer"),
            sdata,
            gene_column="gene",
            ram_threshold_gb=4,
            pixel_size=1 / translation.iloc[0, 0],
        )

    subprocess.run(["rm", "-r", str(zarr_path / "images")])
    subprocess.run(["rm", "-r", str(zarr_path / "points")])

if __name__ == "__main__":
    main()