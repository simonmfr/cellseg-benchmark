import glob
import os
from typing import List, Optional

import geopandas as gpd
import numpy as np
import pandas as pd
import spatialdata as sd
import spatialdata_io
from spatialdata.transformations import Identity, get_transformation, set_transformation

image_based = ["Cellpose", "Negative_Control"]


def process_merscope(sample_name: str, data_dir: str, data_path: str, zmode: str):
    """Load and save a MERSCOPE sample as sdata with specified z_layers configuration. Only loads transcripts and mosaic_images."""
    if zmode not in {"z3", "3d"}:
        raise ValueError(f"Invalid zmode: {zmode}")
    sdata_file = os.path.join(data_dir, "samples", sample_name, f"sdata_{zmode}.zarr")
    if os.path.exists(sdata_file):
        print(f"Skipping {sample_name}: {zmode} file already exists")
        return
    sdata = spatialdata_io.merscope(
        data_path,
        z_layers=3 if zmode == "z3" else range(7),
        backend=None,
        cells_boundaries=False,
        cells_table=False,
        mosaic_images=True,
        transcripts=True,
        slide_name="_".join(sample_name.split("_")[:2]),
        region_name=sample_name.split("_")[2],
    )

    # os.makedirs(os.path.dirname(sdata_file), exist_ok=True) #not necessary
    sdata.write(sdata_file, overwrite=False)
    sdata = sd.read_zarr(sdata_file)

    # set coordinates system
    transformation_to_pixel = get_transformation(
        sdata[list(sdata.points.keys())[0]], "global"
    )

    set_transformation(
        sdata[list(sdata.points.keys())[0]], Identity(), "micron", write_to_sdata=sdata
    )
    set_transformation(
        sdata[list(sdata.points.keys())[0]],
        transformation_to_pixel,
        "pixel",
        write_to_sdata=sdata,
    )

    set_transformation(
        sdata[list(sdata.images.keys())[0]],
        transformation_to_pixel.inverse(),
        "micron",
        write_to_sdata=sdata,
    )
    set_transformation(
        sdata[list(sdata.images.keys())[0]], Identity(), "pixel", write_to_sdata=sdata
    )


def process_merlin_segmentation(
    sample_name, sample_paths, sdata_main, write_to_disk=True
):
    """Process Merlin-specific segmentation data and add it to the main spatial data object.

    Args:
        sample_name: Name of the sample
        sample_paths: Dictionary mapping sample names to their file paths
        sdata_main: Main spatial data object to update
        write_to_disk: Whether to write elements to disk immediately

    Returns:
        Updated sdata_main object
    """
    seg_method = "Cellpose_1_Merlin"

    # Load Merscope data
    if (
        f"boundaries_{seg_method}" not in sdata_main
        or f"adata_{seg_method}" not in sdata_main
    ):
        sdata = spatialdata_io.merscope(
            sample_paths[sample_name],
            z_layers=3,
            backend=None,
            cells_boundaries=True,
            cells_table=True,
            mosaic_images=False,
            transcripts=False,
            slide_name="_".join(sample_name.split("_")[:2]),
            region_name=sample_name.split("_")[2],
        )

    # Handle boundaries
    if f"boundaries_{seg_method}" not in sdata_main:
        polygons = sd.deepcopy(sdata[f"{sample_name}_polygons"])
        sdata_main[f"boundaries_{seg_method}"] = polygons
        if write_to_disk:
            sdata_main.write_element(f"boundaries_{seg_method}")
    else:
        print(f"Skipping {seg_method} as boundaries_{seg_method} exist already.")

    # Handle table
    if f"adata_{seg_method}" not in sdata_main:
        adata = sd.deepcopy(sdata["table"])
        sdata_main[f"adata_{seg_method}"] = adata
        if write_to_disk:
            sdata_main.write_element(f"adata_{seg_method}")
    else:
        print(f"Skipping {seg_method} as adata_{seg_method} exist already.")


def integrate_segmentation_data(
    sdata_path: str,
    seg_methods: List[str],
    sdata_main: sd.SpatialData,
    write_to_disk: bool = True,
    data_path: Optional[str] = None,
):
    """Integrate segmentation data from multiple methods into the main spatial data object.

    Args:
        sdata_path: Path to dircetory of master sdata
        seg_methods: List of segmentation methods to process
        sdata_main: Main spatial data object to update
        write_to_disk: Whether to write elements to disk immediately
        data_path: Optional path to directory to get transformation for adatas

    Returns:
        Updated sdata_main object
    """
    for seg_method in seg_methods:
        seg_path = os.path.join(sdata_path, "results", seg_method, "sdata.zarr")
        if not os.path.exists(os.path.join(seg_path, "shapes")):
            print(f"No boundaries files found for {seg_method}. Skipping.")
            continue
        if not os.path.exists(os.path.join(seg_path, "tables")):
            print(f"No adata files found for {seg_method}. Skipping.")
            continue

        if (
            f"boundaries_{seg_method}" not in sdata_main
            or f"adata_{seg_method}" not in sdata_main
        ):
            sdata = sd.read_zarr(seg_path)

        if f"boundaries_{seg_method}" not in sdata_main:
            # region key specifies the key for boundaries created/used by sopa for this method
            boundary_key = sdata["table"].uns["spatialdata_attrs"]["region"]

            if boundary_key in sdata.shapes.keys():
                sdata_main[f"boundaries_{seg_method}"] = sdata[boundary_key]
                if write_to_disk:
                    sdata_main.write_element(f"boundaries_{seg_method}")
                assign_transformations(sdata_main, seg_method, write_to_disk)
            else:
                print(
                    f"Shapes file missing for {seg_method}. Skipping boundary import. Check conformaty with sopa pipeline, especially sdata['table'].uns['spatialdata_attrs']['region']."
                )
        else:
            print(
                f"Skipping boundary import of {seg_method} as boundaries_{seg_method} exist already."
            )

        # Handle tables
        if f"adata_{seg_method}" not in sdata_main:
            base_dir = os.path.join(seg_path, "tables")
            table_files = [
                os.path.basename(f) for f in glob.glob(os.path.join(base_dir, "table"))
            ]

            if len(table_files) == 1:
                sdata_main[f"adata_{seg_method}"] = sdata[table_files[0]]
                transform_adata(sdata_main, seg_method, data_path=data_path)
                if write_to_disk:
                    sdata_main.write_element(f"adata_{seg_method}")
            elif len(table_files) > 1:
                print(
                    f"Multiple table files found for {seg_method}. Skipping adata import."
                )
            else:
                print(f"Table file missing for {seg_method}. Skipping adata import.")
        else:
            print(
                f"Skipping adata import of {seg_method} as adata_{seg_method} exist already."
            )

    return sdata_main


def assign_transformations(
    sdata_main: sd.SpatialData, seg_method: str, write_to_disk: bool
):
    """Assign transformations to spatial data.

    Args:
        sdata_main: master sdata
        seg_method: current segmentation method
        write_to_disk: if writing to disk
    """
    if write_to_disk:
        backing = sdata_main
    else:
        backing = None

    transformation_to_pixel = get_transformation(
        sdata_main[list(sdata_main.points.keys())[0]], "global"
    )

    if any([seg_method.startswith(method) for method in image_based]):
        set_transformation(
            sdata_main[f"boundaries_{seg_method}"],
            transformation_to_pixel.inverse(),
            "micron",
            write_to_sdata=backing,
        )
        set_transformation(
            sdata_main[f"boundaries_{seg_method}"],
            Identity(),
            "pixel",
            write_to_sdata=backing,
        )
    else:
        set_transformation(
            sdata_main[f"boundaries_{seg_method}"],
            Identity(),
            "micron",
            write_to_sdata=backing,
        )
        set_transformation(
            sdata_main[f"boundaries_{seg_method}"],
            transformation_to_pixel,
            "pixel",
            write_to_sdata=backing,
        )
    return


def transform_adata(sdata_main: sd.SpatialData, seg_method: str, data_path):
    """Add coordinate transforms to adata.

    Args:
        sdata_main: master sdata
        seg_method: current segmentation method
        data_path: path to merscope data
    """
    print(f"Handling adata: {seg_method}")
    transform = pd.read_csv(
        os.path.join(data_path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )

    adata = sdata_main[f"adata_{seg_method}"]
    spatial = adata.obsm["spatial"]

    if any([seg_method.startswith(method) for method in image_based]):
        x = (
            spatial[:, 0] * (1 / transform.iloc[0, 0])
            - (1 / transform.iloc[0, 0]) * transform.iloc[0, 2]
        )
        y = (
            spatial[:, 1] * (1 / transform.iloc[1, 1])
            - (1 / transform.iloc[1, 1]) * transform.iloc[1, 2]
        )
        adata.obsm["spatial_microns"] = np.stack([x, y], axis=1)
        adata.obsm["spatial_pixel"] = spatial
    else:
        adata.obsm["spatial_microns"] = spatial
        x = spatial[:, 0] * transform.iloc[0, 0] + transform.iloc[0, 2]
        y = spatial[:, 1] * transform.iloc[1, 1] + transform.iloc[1, 2]
        adata.obsm["spatial_pixel"] = np.stack([x, y], axis=1)
    adata.uns["spatialdata_attrs"]["region"] = f"boundaries_{seg_method}"
    adata.uns["spatialdata_attrs"]["region_key"] = "region"
    adata.obs["region"] = f"boundaries_{seg_method}"
    adata.obs["region"] = pd.Categorical(adata.obs["region"])
    adata.obs["cell_id"] = adata.obs[adata.uns["spatialdata_attrs"]["instance_key"]]
    adata.obs["cell_id"] = adata.obs["cell_id"]
    adata.uns["spatialdata_attrs"]["instance_key"] = "cell_id"

    sdata_main[f"adata_{seg_method}"] = adata
    return


def update_element(sdata, element_name):
    """Workaround for updating a backed element in sdata.

    Adapted from https://github.com/scverse/spatialdata/blob/main/tests/io/test_readwrite.py#L156
    """
    new_name = f"{element_name}_tmp"
    name = element_name
    # a a. write a backup copy of the data
    sdata[new_name] = sdata[name]
    sdata.write_element(new_name)
    # a2. remove the in-memory copy from the SpatialData object (note,
    # at this point the backup copy still exists on-disk)
    del sdata[new_name]
    del sdata[name]
    # a3 load the backup copy into memory
    sdata_copy = sd.read_zarr(sdata.path)
    # b1. rewrite the original data
    sdata.delete_element_from_disk(name)
    sdata[name] = sdata_copy[new_name]
    sdata.write_element(name, overwrite=True)
    # b2. reload the new data into memory (because it has been written but in-memory it still points
    # from the backup location)
    sdata = sd.read_zarr(sdata.path)
    # c. remove the backup copy
    del sdata[new_name]
    sdata.delete_element_from_disk(new_name)


def pixel_to_microns(
    sdata: sd.SpatialData,
    transform_file: str,
    shape_patterns: List[str] = None,
    exclude_patterns: List[str] = None,
    overwrite: bool = False,
):
    """Transform Cellpose boundaries from pixels to microns.

    Parameters:
    -----------
    sdata : SpatialData
        SpatialData object containing shape collections
    transform_file : str
        Path to file containing transformation matrix
    shape_patterns : list, optional
        List of patterns to match for shape names
    exclude_patterns : list, optional
        List of patterns to exclude from shape names
    overwrite : bool, optional
        Whether to overwrite existing shapes (default: False)
    """
    if exclude_patterns is None:
        exclude_patterns = []
    if shape_patterns is None:
        shape_patterns = sdata.shapes.keys()

    transform_matrix = pd.read_csv(transform_file, sep=" ", header=None).values
    inv_matrix = np.linalg.inv(transform_matrix)

    # Track number of transformations
    transform_count = 0

    # Process eligible shapes
    for shape_name in list(sdata.shapes.keys()):
        # Check if shape name matches any of the patterns and doesn't match any exclude patterns
        if any(pattern in shape_name for pattern in shape_patterns) and not any(
            exclude in shape_name for exclude in exclude_patterns
        ):
            # Get boundaries
            boundaries = sdata.shapes[shape_name]

            # Create transformed geometries
            transformed_geoms = boundaries.affine_transform(
                [
                    inv_matrix[0, 0],
                    inv_matrix[0, 1],
                    inv_matrix[1, 0],
                    inv_matrix[1, 1],
                    inv_matrix[0, 2],
                    inv_matrix[1, 2],
                ]
            )

            # Create new GeoDataFrame
            transformed_boundaries = gpd.GeoDataFrame(
                boundaries.drop(columns=["geometry"]),
                geometry=transformed_geoms,
                crs=boundaries.crs,
            )

            # Copy attributes
            for key, value in boundaries.attrs.items():
                transformed_boundaries.attrs[key] = value

            # Determine output name
            output_name = shape_name if overwrite else f"{shape_name}_microns"

            # Add to spatialdata
            sdata.shapes[output_name] = transformed_boundaries
            transform_count += 1
