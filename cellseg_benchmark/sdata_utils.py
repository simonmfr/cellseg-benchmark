import gzip
import io
import logging
import os
from os import listdir
from os.path import join
from re import split
from typing import List, Optional, Union, Dict

import geopandas as gpd
import numpy as np
import pandas as pd
import spatialdata as sd
import spatialdata_io
from spatialdata.transformations import Identity, get_transformation, set_transformation
from tifffile import imread
from tqdm import tqdm

from ._constants import image_based, methods_3D
from .ficture_utils import create_factor_level_image, parse_metadata


def process_merscope(sample_name: str, data_dir: str, data_path: str, zmode: str) -> None:
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
    sample_name: str, sample_paths: Dict[str, str], sdata_main: sd.SpatialData, write_to_disk: bool=True
)-> None:
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
    logger: Optional[logging.Logger] = None,
) -> sd.SpatialData:
    """Integrate segmentation data from multiple methods into the main spatial data object.

    Args:
        sdata_path: Path to dircetory of master sdata
        seg_methods: List of segmentation methods to process
        sdata_main: Main spatial data object to update
        write_to_disk: Whether to write elements to disk immediately
        data_path: Optional path to directory to get transformation for adatas
        logger: Optional logger object to write messages to console

    Returns:
        Updated sdata_main object
    """
    for seg_method in tqdm(seg_methods):
        if logger is not None:
            logger.info(f"Adding {seg_method}...")
        seg_path = os.path.join(sdata_path, "results", seg_method, "sdata.zarr")
        if not os.path.exists(os.path.join(seg_path, "shapes")):
            if logger:
                logger.warning(
                    "No boundaries files found for {}. Skipping".format(seg_method)
                )
            else:
                print(f"No boundaries files found for {seg_method}. Skipping.")
            continue
        if not os.path.exists(os.path.join(seg_path, "tables")):
            if logger:
                logger.warning(
                    "No adata files found for {}. Skipping".format(seg_method)
                )
            else:
                print(f"No adata files found for {seg_method}. Skipping.")
            continue

        if (
            f"boundaries_{seg_method}" not in sdata_main
            or f"adata_{seg_method}" not in sdata_main
        ):
            sdata = sd.read_zarr(seg_path)

        if f"boundaries_{seg_method}" not in sdata_main:
            sdata_main = build_shapes(
                sdata, sdata_main, seg_method, sdata_path, write_to_disk=write_to_disk
            )
        else:
            if logger:
                logger.warning(
                    "Skipping boundary import of {} as boundaries_{} exist already.".format(
                        seg_method, seg_method
                    )
                )
            else:
                print(
                    f"Skipping boundary import of {seg_method} as boundaries_{seg_method} exist already."
                )

        # Handle tables
        if f"adata_{seg_method}" not in sdata_main:
            if len(sdata.tables) == 1:
                sdata_main[f"adata_{seg_method}"] = sdata[
                    list(sdata.tables.keys())[0]
                ].copy()
                transform_adata(sdata_main, seg_method, data_path=data_path)

                if os.path.exists(
                    join(
                        sdata_path,
                        "results",
                        seg_method,
                        "cell_type_annotation",
                        "adata_obs_annotated.csv",
                    )
                ):
                    sdata_main = add_cell_type_annotation(
                        sdata_main, sdata_path, seg_method, write_to_disk=write_to_disk
                    )
                else:
                    if logger:
                        logger.warning(
                            "No annotation files found for {}. Skipping".format(
                                seg_method
                            )
                        )
                    else:
                        print(
                            f"No annotation files found for {seg_method}. Skipping annotation."
                        )
                if "volume" not in sdata_main[f"adata_{seg_method}"].obs.columns:
                    sdata_main = calculate_volume(
                        seg_method, sdata_main, sdata_path, write_to_disk=write_to_disk
                    )
                if os.path.exists(
                    join(sdata_path, "results", seg_method, "Ficture_stats")
                ):
                    logger.info("Adding Ficture stats to {}...".format(seg_method))
                    add_ficture(
                        sdata_main, seg_method, sdata_path, write_to_disk=write_to_disk
                    )
                else:
                    logger.warning(
                        "No Ficture_stats files found for {}. Skipping.".format(
                            seg_method
                        )
                    )
            elif len(sdata.tables) > 1:
                if logger:
                    logger.warning(
                        "Multiple adata files found for {}. Skipping adata import".format(
                            seg_method
                        )
                    )
                else:
                    print(
                        f"Multiple table files found for {seg_method}. Skipping adata import."
                    )
            else:
                if logger:
                    logger.warning(
                        "No adata files found for {}. Skipping adata import".format(
                            seg_method
                        )
                    )
                else:
                    print(
                        f"Table file missing for {seg_method}. Skipping adata import."
                    )
        else:
            if logger:
                logger.warning(
                    "Skipping adata import of {} as adata_{} exist already.".format(
                        seg_method, seg_method
                    )
                )
            else:
                print(
                    f"Skipping adata import of {seg_method} as adata_{seg_method} exist already."
                )

    return sdata_main


def build_shapes(sdata: sd.SpatialData, sdata_main: sd.SpatialData, seg_method: str, sdata_path: str,
                 write_to_disk: bool, logger: logging.Logger=None):
    """Insert shapes of segmentation method into sdata_main."""
    boundary_key = sdata["table"].uns["spatialdata_attrs"]["region"]

    if seg_method.startswith("Proseg"):
        path = join(
            sdata_path,
            "results",
            seg_method,
            "sdata.zarr",
            ".sopa_cache",
            "transcript_patches",
            "0",
            "cell-polygons-layers.geojson.gz",
        )
        with gzip.open(path, "rt", encoding="utf-8") as f:
            geojson_text = f.read()
        geojson_io = io.StringIO(geojson_text)
        gdf = gpd.read_file(geojson_io)
        gdf = gdf.merge(sdata["table"].obs[["cell", "cell_id"]], on="cell")
        sdata_main[f"boundaries_{seg_method}"] = gdf
    elif boundary_key in sdata.shapes.keys():
        sdata_main[f"boundaries_{seg_method}"] = sdata[boundary_key]
        if write_to_disk:
            if (
                join("shapes", f"boundaries_{seg_method}")
                in sdata_main.elements_paths_on_disk()
            ):
                update_element(sdata_main, f"boundaries_{seg_method}")
            else:
                sdata_main.write_element(f"boundaries_{seg_method}")
        assign_transformations(sdata_main, seg_method, write_to_disk)
    else:
        if logger:
            logger.warning(
                "Shapes file missing for {}. Skipping boundary import. Check conformaty with sopa pipeline, especially sdata['table'].uns['spatialdata_attrs']['region'].".format(
                    seg_method
                )
            )
        else:
            print(
                f"Shapes file missing for {seg_method}. Skipping boundary import. Check conformaty with sopa pipeline, especially sdata['table'].uns['spatialdata_attrs']['region']."
            )
    return sdata_main


def add_cell_type_annotation(
    sdata_main: sd.SpatialData, sdata_path: str, seg_method: str, write_to_disk: bool, logger: logging.Logger=None
) -> sd.SpatialData:
    """Add cell type annotations to sdata_main, including adding volumes."""
    cell_type_information = [
        "cell_type_mmc_incl_low_quality_revised",
        "cell_type_mmc_incl_low_quality_clusters",
        "cell_type_mmc_incl_low_quality",
        "cell_type_mmc_incl_mixed_revised",
        "cell_type_mmc_incl_mixed_clusters",
        "cell_type_mmc_incl_mixed",
        "cell_type_mmc_raw_revised",
        "cell_type_mmc_raw_clusters",
        "cell_type_mmc_raw",
        "cell_id",
    ]
    try:
        cell_type = pd.read_csv(
            join(
                sdata_path,
                "results",
                seg_method,
                "cell_type_annotation",
                "adata_obs_annotated.csv",
            )
        )[cell_type_information]
    except KeyError:
        if logger:
            logger.warning(
                "no propper cell annotation found for {}. Skipping.".format(seg_method)
            )
        return sdata_main
    if set(cell_type_information) & set(sdata_main[f"adata_{seg_method}"].obs.columns):
        for col in set(cell_type_information) & set(
            sdata_main[f"adata_{seg_method}"].obs.columns
        ):
            del sdata_main[f"adata_{seg_method}"].obs[col]
    new_obs = sdata_main[f"adata_{seg_method}"].obs.merge(
        cell_type, how="left", left_index=True, right_on="cell_id"
    )
    new_obs.index = sdata_main[f"adata_{seg_method}"].obs.index
    for col in new_obs.columns:
        if isinstance(new_obs[col].dtype, pd.CategoricalDtype):
            new_obs[col] = new_obs[col].cat.add_categories("Low-Read-Cells")
        new_obs[col].fillna("Low-Read-Cells", inplace=True)
    sdata_main[f"adata_{seg_method}"].obs = new_obs
    if write_to_disk:
        if join("tables", f"adata_{seg_method}") in sdata_main.elements_paths_on_disk():
            update_element(sdata_main, f"adata_{seg_method}")
        else:
            sdata_main.write_element(f"adata_{seg_method}")
    return sdata_main


def add_ficture(sdata_main: sd.SpatialData, seg_method: str, sdata_path: str, write_to_disk: bool) -> sd.SpatialData:
    """Add ficture information to sdata_main."""
    adata = sdata_main[f"adata_{seg_method}"]
    for file in os.listdir(join(sdata_path, "results", seg_method, "Ficture_stats")):
        name = file.split(".")[0]
        adata.obsm[name] = pd.read_csv(
            join(sdata_path, "results", seg_method, "Ficture_stats", file), index_col=0
        )
    sdata_main[f"adata_{seg_method}"] = adata
    if write_to_disk:
        if join("tables", f"adata_{seg_method}") in sdata_main.elements_paths_on_disk():
            update_element(sdata_main, f"adata_{seg_method}")
        else:
            sdata_main.write_element(f"adata_{seg_method}")
    return sdata_main


def calculate_volume(
    seg_method: str, sdata_main: sd.SpatialData, sdata_path: str, write_to_disk=False, logger=None
) -> sd.SpatialData:  #
    """Calculate volume of sdata."""
    adata = sdata_main[f"adata_{seg_method}"]
    if any([seg_method.startswith(x) for x in methods_3D]):
        if seg_method.startswith("Proseg"):  # boundaries in microns
            path = join(
                sdata_path,
                "results",
                seg_method,
                "sdata.zarr",
                ".sopa_cache",
                "transcript_patches",
                "0",
                "cell-polygons-layers.geojson.gz",
            )
            with gzip.open(path, "rt", encoding="utf-8") as f:
                geojson_text = f.read()

            geojson_io = io.StringIO(geojson_text)

            gdf = gpd.read_file(geojson_io)
            counts = gdf.cell
            ind = list(gdf.cell)
            counts.index = ind
            counts = counts.groupby(level=0).count()
            slice_height = (
                10 / counts.max()
            )  # 10 um total. Assume at least one detected cell stretches whole slice
            gdf["volume_per_slice"] = [x.area * slice_height for x in gdf["geometry"]]
            area = gdf[["cell", "volume_per_slice"]].groupby("cell").sum()
            area.rename(columns={"volume_per_slice": "volume"}, inplace=True)
            if "volume" in adata.obs.columns:
                del adata.obs["volume"]
            tmp = adata.obs.merge(area, how="left", left_on="cell", right_on="cell")
            tmp.index = adata.obs.index
            adata.obs = tmp
        elif seg_method.startswith("vpt_3D"):
            boundaries = sdata_main.transform_element_to_coordinate_system(
                f"boundaries_{seg_method}", target_coordinate_system="micron"
            )
            slice_height = 10 / boundaries[["ZIndex"]].groupby(level=0).count().max()
            boundaries["volume_per_slice"] = [
                x.area * slice_height for x in boundaries["Geometry"]
            ]
            area = (
                boundaries[["EntityID", "volume_per_slice"]].groupby("EntityID").sum()
            )
            area.rename(columns={"volume_per_slice": "volume"}, inplace=True)
            if "volume" in adata.obs.columns:
                del adata.obs["volume"]
            tmp = adata.obs.merge(
                area, how="left", left_index=True, right_on="EntityID"
            )
            tmp.index = adata.obs.index
            adata.obs = tmp
    else:
        try:
            boundaries = sdata_main.transform_element_to_coordinate_system(
                f"boundaries_{seg_method}", target_coordinate_system="micron"
            )
        except KeyError:
            if logger:
                logger.warning(
                    "Volume cannot be computed for {}. Skipping.".format(seg_method)
                )
            return sdata_main
        adata.obs["volume"] = boundaries.geometry.area * 10
    sdata_main[f"adata_{seg_method}"] = adata
    if write_to_disk:
        if join("tables", f"adata_{seg_method}") in sdata_main.elements_paths_on_disk():
            update_element(sdata_main, f"adata_{seg_method}")
        else:
            sdata_main.write_element(f"adata_{seg_method}")
    return sdata_main


def assign_transformations(
    sdata_main: sd.SpatialData, seg_method: str, write_to_disk: bool
) -> None:
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


def transform_adata(sdata_main: sd.SpatialData, seg_method: str, data_path: str) -> None:
    """Add coordinate transforms to adata.

    Args:
        sdata_main: master sdata
        seg_method: current segmentation method
        data_path: path to merscope data
    """
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


def update_element(sdata: sd.SpatialData, element_name: str) -> None:
    """Workaround for updating a backed element in sdata.

    Adapted from https://github.com/scverse/spatialdata/blob/main/tests/io/test_readwrite.py#L156
    """
    new_name = f"{element_name}_tmp"
    name = element_name
    # a a. write a backup copy of the data
    sdata[new_name] = sdata[name]
    sdata.write_element(new_name, overwrite=True)
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
) -> None:
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


def prepare_ficture(
    data_path: str,
    results_path: str,
    top_n_factors: int=3,
    n_ficture: int=21,
    logger: logging.Logger=None,
    factors: Optional[List[int]] = None,
) -> Dict[str, Union[np.ndarray, List[int]]]:
    """Generate ficture images stack and other ficture information.

    Args:
        data_path: Path to merscope data
        results_path: path to Segmentation folder
        n_ficture: number of factors of ficture run
        logger: logger instance
        factors: if provided, only these ficture images will be generated.

    Returns:
        ficture images and factors of the images
    """
    if logger is not None:
        logger.info(f"Generating ficture images for {data_path}")
    DAPI_shape = imread(join(data_path, "images/mosaic_DAPI_z3.tif")).shape
    transform = pd.read_csv(
        join(data_path, "images/micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
    )

    if "Ficture" not in listdir(results_path):
        return {}
    ficture_path = join(results_path, "Ficture", "output")
    for file in listdir(ficture_path):
        if n_ficture == int(split(r"\.|F", file)[1]):
            ficture_path = join(ficture_path, file)

    ficture_full_path = ""
    for file in listdir(ficture_path):
        if file.endswith(".pixel.sorted.tsv.gz"):
            ficture_full_path = join(ficture_path, file)
    assert ficture_full_path != "", "Ficture output not correctly computed."

    fic_header = ["BLOCK", "X", "Y", "K1", "K2", "K3", "P1", "P2", "P3"]
    ficture_pixels = pd.read_csv(
        ficture_full_path, sep="\t", names=fic_header, comment="#"
    )

    metadata = parse_metadata(ficture_full_path)
    scale = float(metadata["SCALE"])
    offset_x = float(metadata["OFFSET_X"])
    offset_y = float(metadata["OFFSET_Y"])
    #assume, that transform[0,1], transform[1,0] = 0
    ficture_pixels["X_pixel"] = (
        ficture_pixels["X"] / scale * transform.iloc[0, 0]
        + offset_x * transform.iloc[0, 0]
        + transform.iloc[0, 2]
    )
    ficture_pixels["Y_pixel"] = (
        ficture_pixels["Y"] / scale * transform.iloc[1, 1]
        + offset_y * transform.iloc[1, 1]
        + transform.iloc[1, 2]
    )
    del transform, metadata

    unique_factors = set()
    for i in range(1, top_n_factors + 1):
        unique_factors= unique_factors.union(set(np.unique(ficture_pixels[f"K{i}"])))
    unique_factors = list(unique_factors)

    if factors is not None:
        assert all([x in unique_factors for x in factors])
        unique_factors = list(set(factors))

    for factor in tqdm(unique_factors):
        if logger is not None:
            logger.info(f"Building ficture image for {factor}")
        try:
            image_stack
        except NameError:
            image_stack = create_factor_level_image(ficture_pixels, factor, DAPI_shape, top_n_factors)
        else:
            image_stack = np.concatenate(
                (
                    image_stack,
                    create_factor_level_image(ficture_pixels, factor, DAPI_shape, top_n_factors),
                ),
                axis=0,
                dtype=np.uint16,
            )
    return {"images":image_stack, "factors": unique_factors}
