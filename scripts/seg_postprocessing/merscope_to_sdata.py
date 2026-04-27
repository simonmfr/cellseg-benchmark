import argparse
import pathlib
import subprocess
import pandas
import spatialdata
import spatialdata_io
import sopa.aggregation
import sopa.io.explorer

def main(args):
    """Convert default Merscope segmentation output to sdata."""
    data_path = pathlib.Path(args.data_path)
    save_path = pathlib.Path(args.save_path)

    sdata = spatialdata_io.merscope(
        data_path,
        vpt_outputs={
            "cell_by_gene": data_path / "cell_by_gene.csv",
            "cell_metadata": data_path / "cell_metadata.csv",
            "cell_boundaries": data_path / "cell_boundaries.parquet",
        },
    )
    sdata["table"].obs.rename(columns={"EntityID": "cell_id"}, inplace=True)
    shapes_key = list(sdata.shapes.keys())[0]
    boundaries = sdata[shapes_key]
    boundaries.rename(columns={"EntityID": "cell_id"}, inplace=True)
    boundaries.index = boundaries["cell_id"]
    sdata["table"].uns["spatialdata_attrs"]["instance_key"] = "cell_id"
    table_ids = sdata["table"].obs["cell_id"].astype(str)
    boundary_ids = boundaries.index.astype(str)
    assert len(table_ids) == len(boundary_ids), "Number of table and boundary ids do not match."
    if not (table_ids.values == boundary_ids.values).all():
        boundaries = boundaries.loc[table_ids]
        sdata[shapes_key] = boundaries
    translation = pandas.read_csv(
        data_path / "images" / "micron_to_mosaic_pixel_transform.csv",
        sep=" ",
        header=None,
    )
    sdata["table"].obsm["intensities"] = sopa.aggregation.aggregate_channels(sdata, shapes_key=shapes_key)
    sdata.write(str(save_path / "sdata.zarr"), overwrite=True)
    sdata = spatialdata.read_zarr(str(save_path / "sdata.zarr"))
    if args.explorer:
        sopa.io.explorer.write(
            str(save_path / "sdata.explorer"),
            sdata,
            gene_column="gene",
            ram_threshold_gb=4,
            pixel_size=1/translation.iloc[0, 0],
        )
    subprocess.run(["rm", "-r", str(save_path / "sdata.zarr" / "images")])
    subprocess.run(["rm", "-r", str(save_path / "sdata.zarr" / "points")])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert default Merscope segmentation output to sdata."
    )
    parser.add_argument("data_path", help="Path to data folder.")
    parser.add_argument("save_path", help="Path to output folder.")
    parser.add_argument(
        "--explorer", action="store_true", help="if to compute 10X Xenium explorer files"
    )
    args = parser.parse_args()
    main(args)