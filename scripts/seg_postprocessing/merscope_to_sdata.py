import argparse
import pathlib
import subprocess
import pandas
import geopandas
import shapely
import spatialdata
import spatialdata_io
import spatialdata.models
import sopa.aggregation
import sopa.io.explorer

def main():
    parser = argparse.ArgumentParser(description="Convert Merscope segmentation output to sdata.")
    parser.add_argument("data_path", help="Path to merscope folder.")
    parser.add_argument("save_path", help="Path to output folder.")
    parser.add_argument(
        "--segmentation",
        choices=["cellpose", "watershed"],
        default="cellpose",
        help="Segmentation type. cellpose=2D, watershed=pseudo-3D (one polygon per z-plane).",
    )
    parser.add_argument(
        "--boundaries-parquet",
        default=None,
        help="Path to boundaries parquet. Defaults to <data_path>/cell_boundaries.parquet.",
    )
    parser.add_argument("--explorer", action="store_true", help="Whether to compute 10X Xenium explorer file (cellpose only).")
    args = parser.parse_args()

    if args.explorer and args.segmentation != "cellpose":
        parser.error("--explorer is only supported for cellpose (2D) segmentation.")

    data_path = pathlib.Path(args.data_path)
    save_path = pathlib.Path(args.save_path)
    boundaries_path = pathlib.Path(args.boundaries_parquet) if args.boundaries_parquet else data_path / "cell_boundaries.parquet"

    sdata = spatialdata_io.merscope(
        data_path,
        vpt_outputs={
            "cell_by_gene": data_path / "cell_by_gene.csv",
            "cell_metadata": data_path / "cell_metadata.csv",
            "cell_boundaries": boundaries_path,
        },
    )

    sdata["table"].obs.rename(columns={"EntityID": "cell_id"}, inplace=True)
    sdata["table"].uns["spatialdata_attrs"]["instance_key"] = "cell_id"
    table_ids = sdata["table"].obs["cell_id"].astype(str)

    if args.segmentation == "cellpose":
        shapes_key = list(sdata.shapes.keys())[0]
        boundaries = sdata[shapes_key]
        boundaries.rename(columns={"EntityID": "cell_id"}, inplace=True)
        boundaries.index = boundaries["cell_id"].astype(str)
        boundaries = boundaries.loc[table_ids.values]
        sdata[shapes_key] = boundaries
        agg_key = shapes_key
    elif args.segmentation == "watershed":
        boundaries = geopandas.read_parquet(
            boundaries_path,
            columns=("ID", "EntityID", "ZIndex", "ZLevel", "Geometry"),
        )
        boundaries.rename_geometry("geometry", inplace=True)
        boundaries.rename(columns={"EntityID": "cell_id"}, inplace=True)
        boundaries.set_index("cell_id", drop=False, inplace=True)
        boundaries.index = boundaries.index.rename(None)

        boundaries["geometry"] = boundaries["geometry"].make_valid()
        n_collections = (boundaries["geometry"].geom_type == "GeometryCollection").sum()
        if n_collections > 0:
            print(f"Warning: {n_collections} cells had invalid geometries after make_valid(), polygon parts will be extracted.")
        boundaries["geometry"] = boundaries["geometry"].apply(
            lambda g: shapely.geometry.MultiPolygon([p for p in g.geoms if p.geom_type == "Polygon"])
            if g.geom_type == "GeometryCollection" else g
        )
        boundaries_2d = boundaries[["cell_id", "geometry"]].dissolve(by="cell_id")
        boundaries_2d.index = boundaries_2d.index.rename(None)

        sdata["boundaries_vpt_3D"] = spatialdata.models.ShapesModel.parse(boundaries)
        sdata["boundaries_vpt_2D"] = spatialdata.models.ShapesModel.parse(boundaries_2d)

        sdata["table"].uns["spatialdata_attrs"]["region"] = "boundaries_vpt_3D"
        region_key = sdata["table"].uns["spatialdata_attrs"]["region_key"]
        sdata["table"].obs[region_key] = "boundaries_vpt_3D"
        sdata["table"].obs[region_key] = sdata["table"].obs[region_key].astype("category")

        agg_key = "boundaries_vpt_2D"

    translation = pandas.read_csv(
        data_path / "images" / "micron_to_mosaic_pixel_transform.csv",
        sep=" ",
        header=None,
    )

    sdata["table"].obsm["intensities"] = sopa.aggregation.aggregate_channels(sdata, shapes_key=agg_key)

    if args.segmentation == "watershed":
        del sdata["boundaries_vpt_2D"]

    sdata.write(str(save_path / "sdata.zarr"), overwrite=True)
    sdata = spatialdata.read_zarr(str(save_path / "sdata.zarr"))

    if args.explorer:
        sopa.io.explorer.write(
            str(save_path / "sdata.explorer"),
            sdata,
            gene_column="gene",
            ram_threshold_gb=4,
            pixel_size=1 / translation.iloc[0, 0],
        )

    subprocess.run(["rm", "-r", str(save_path / "sdata.zarr" / "images")])
    subprocess.run(["rm", "-r", str(save_path / "sdata.zarr" / "points")])

if __name__ == "__main__":
    main()