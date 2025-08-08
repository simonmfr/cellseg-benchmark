import ast
import gzip
import io
import logging
import math
import os
import warnings
from os import listdir
from os.path import join
from re import split
from typing import Dict, List, Optional, Union

import geopandas as gpd
import numpy as np
import pandas as pd
import spatialdata as sd
import spatialdata_io
from joblib import Parallel, delayed
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon
from shapely.ops import unary_union
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity, get_transformation, set_transformation, Affine
from tifffile import imread
from tqdm import tqdm

from ._constants import image_based, methods_3D
from .ficture_utils import create_factor_level_image, parse_metadata

PI = math.pi

def process_merscope(
    sample_name: str, data_dir: str, data_path: str, zmode: str
) -> None:
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
    sample_name: str,
    sample_paths: Dict[str, str],
    sdata_main: sd.SpatialData,
    write_to_disk: bool = True,
) -> None:
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
                sdata,
                sdata_main,
                seg_method,
                sdata_path,
                logger=logger,
            )
            if write_to_disk:
                sdata_main.write_element(f"boundaries_{seg_method}")
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
            if logger:
                logger.info(f"Adding adata for {seg_method}...")
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
                        sdata_main,
                        sdata_path,
                        seg_method,
                        logger=logger,
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
                if "volume_final" not in sdata_main[f"adata_{seg_method}"].obs.columns:
                    if any([seg_method.startswith(x) for x in methods_3D]):
                        n_planes_2d = None
                    else:
                        n_planes_2d = 7
                    sdata_main = calculate_volume(
                        seg_method,
                        sdata_main,
                        write_to_disk=write_to_disk,
                        n_planes_2d=n_planes_2d,
                        logger=logger,
                    )
                if os.path.exists(
                    join(sdata_path, "results", seg_method, "Ficture_stats")
                ):
                    logger.info("Adding Ficture stats to {}...".format(seg_method))
                    add_statistical_data(
                        sdata_main, seg_method, sdata_path
                    )
                else:
                    logger.warning(
                        "No Ficture_stats files found for {}. Skipping.".format(
                            seg_method
                        )
                    )
                if write_to_disk:
                    sdata_main.write_element(f"adata_{seg_method}")
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


def build_shapes(
    sdata: sd.SpatialData,
    sdata_main: sd.SpatialData,
    seg_method: str,
    sdata_path: str,
    logger: logging.Logger = None,
):
    """Insert shapes of segmentation method into sdata_main."""
    if logger:
        logger.info(f"Adding shapes of {seg_method}...")
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
        sdata_main[f"boundaries_{seg_method}"] = ShapesModel.parse(gdf)
        assign_transformations(sdata_main, seg_method)
    elif boundary_key in sdata.shapes.keys():
        sdata_main[f"boundaries_{seg_method}"] = sdata[boundary_key]
        assign_transformations(sdata_main, seg_method)
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
    sdata_main: sd.SpatialData,
    sdata_path: str,
    seg_method: str,
    logger: logging.Logger = None,
) -> sd.SpatialData:
    """Add cell type annotations to sdata_main, including adding volumes."""
    if logger:
        logger.info(f"Adding cell type annotations for {seg_method}...")
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
    try:
        if "cell_id" in sdata_main[f"adata_{seg_method}"].obs.columns:
            sdata_main[f"adata_{seg_method}"].obs.drop(columns="cell_id", inplace=True)
        new_obs = sdata_main[f"adata_{seg_method}"].obs.merge(
            cell_type, how="left", left_index=True, right_on="cell_id"
        )
    except ValueError:
        cell_type["cell_id"] = cell_type["cell_id"].astype(str)
        new_obs = sdata_main[f"adata_{seg_method}"].obs.merge(
            cell_type, how="left", left_index=True, right_on="cell_id"
        )
    new_obs.index = sdata_main[f"adata_{seg_method}"].obs.index
    for col in new_obs.columns:
        if isinstance(new_obs[col].dtype, pd.CategoricalDtype):
            new_obs[col] = new_obs[col].cat.add_categories("Low-Read-Cells")
        new_obs[col].fillna("Low-Read-Cells", inplace=True)
    sdata_main[f"adata_{seg_method}"].obs = new_obs
    return sdata_main


def add_statistical_data(
    sdata_main: sd.SpatialData, seg_method: str, sdata_path: str
) -> sd.SpatialData:
    """Add ficture and ovrlpy information to sdata_main."""
    adata = sdata_main[f"adata_{seg_method}"]
    for file in os.listdir(join(sdata_path, "results", seg_method, "Ficture_stats")):
        name = file.split(".")[0]
        ficture_stats = pd.read_csv(
            join(sdata_path, "results", seg_method, "Ficture_stats", file), index_col=0
        )
        ficture_stats.index = ficture_stats.index.astype(str)
        adata.obsm[f"ficture_{name}"] = ficture_stats
    for file in os.listdir(join(sdata_path, "results", seg_method, "Ovrlpy_stats")):
        if file.endswith(".csv"):
            name = file.split(".")[0]
            ovrlpy_stats = pd.read_csv(
                join(sdata_path, "results", seg_method, "Ovrlpy_stats", file), index_col=0
            )
            ovrlpy_stats.index = ovrlpy_stats.index.astype(str)
            adata.obsm[name] = ovrlpy_stats
    sdata_main[f"adata_{seg_method}"] = adata
    return sdata_main


def calculate_volume(
    seg_method: str,
    sdata_main: sd.SpatialData,
    z_spacing: float = 1.5,
    n_planes_2d: Optional[int] = None,
    logger: logging.Logger = None,
):
    """Calculate volume of cells in sdata_main."""
    boundaries = sd.transform(
        sdata_main[f"boundaries_{seg_method}"], to_coordinate_system="micron"
    )
    if any([seg_method.startswith(x) for x in methods_3D]):
        if seg_method.startswith("Proseg"):
            z_level_name = "layer"
            cell_identifier = "cell_id"
        elif seg_method.startswith("vpt_3D"):
            z_level_name = "ZIndex"
            cell_identifier = "cell_id"
        if logger:
            logger.info(f"collecting volume metadata for {seg_method}")
        global_z_min, global_z_max = (
            boundaries[z_level_name].min(),
            boundaries[z_level_name].max(),
        )
        try:
            grouped = boundaries.groupby(cell_identifier)
        except ValueError:
            boundaries.index.rename(None, inplace=True)
            grouped = boundaries.groupby(cell_identifier)

        items = [
            (entity_id, group[[z_level_name, "geometry"]])
            for entity_id, group in grouped
        ]
        morphology_rows = Parallel(n_jobs=-1, backend="loky")(
            delayed(compute_polygon_stats_3D)(eid, grp, z_spacing, z_level_name, global_z_min, global_z_max, logger) for eid, grp in items
        )

    else:
        assert n_planes_2d is not None, "provide n_planes_2d parameter for 2D methods"
        grouped = boundaries.groupby(level=0)

        if logger:
            logger.info(f"calculate volume metrics {seg_method}")
        scale = z_spacing * n_planes_2d
        items = [
            (entity_id, group.geometry.iat[0])
            for entity_id, group in grouped
        ]
        morphology_rows = Parallel(n_jobs=-1, backend="loky")(
            delayed(compute_polygon_stats_2D)(item, scale, logger) for item in items
        )
        morphology_rows = [r for r in morphology_rows if r is not None]

    if morphology_rows:
        df_morph = pd.DataFrame(morphology_rows)
        df_morph["cell_id"] = df_morph["cell_id"].astype(str)
        df_morph.index = df_morph["cell_id"]
        if "dimensionality" in df_morph.columns:
            df_morph.drop(columns=["dimensionality", "cell_id"], inplace=True)
        else:
            df_morph.drop(columns=["cell_id"], inplace=True)

        adata = sdata_main.tables[f"adata_{seg_method}"]
        if "area" in adata.obs.columns:
            adata.obs.drop(columns="area", inplace=True)
        if "solidity" in adata.obs.columns:
            adata.obs.drop(columns="solidity", inplace=True)
        adata.obs = adata.obs.merge(df_morph, left_index=True, right_index=True)
        sdata_main[f"adata_{seg_method}"] = adata
    return sdata_main


def compute_polygon_stats_2D(entity_geom_pair, scale, logger=None):
    entity_id, geom = entity_geom_pair
    try:
        data = _compute_2d_metrics(geom, scale)
        data["cell_id"] = entity_id
        return data
    except Exception as e:
        if logger is not None:
            logger.warning(f"Failed to process entity {entity_id}: {str(e)}")
        else:
            print(f"Failed {entity_id}: {e}")
        return None

def compute_polygon_stats_3D(entity_id, group, z_spacing, z_level_name, global_z_min, global_z_max, logger=None):
    try:
        m = _compute_3d_metrics(
            group,
            z_spacing,
            global_z_min=global_z_min,
            global_z_max=global_z_max,
            ZIndex=z_level_name,
        )
        m["cell_id"] = entity_id
        return m
    except Exception as e:
        if logger is not None:
            logger.warning(f"Failed to process entity {entity_id}: {str(e)}")
        else:
            print(f"Failed {entity_id}: {e}")
        return None


def assign_transformations(
    sdata_main: sd.SpatialData, seg_method: str
) -> None:
    """Assign transformations to spatial data.

    Args:
        sdata_main: master sdata
        seg_method: current segmentation method
        write_to_disk: if writing to disk
    """
    transformation_to_pixel = get_transformation(
        sdata_main[list(sdata_main.points.keys())[0]], "global"
    )

    if any([seg_method.startswith(method) for method in image_based]):
        if seg_method == "Cellpose_1_Merlin":
            set_transformation(
                sdata_main[f"boundaries_{seg_method}"],
                Identity(),
                "micron"
            )
            set_transformation(
                sdata_main[f"boundaries_{seg_method}"],
                transformation_to_pixel,
                "pixel"
            )
        else:
            set_transformation(
                sdata_main[f"boundaries_{seg_method}"],
                transformation_to_pixel.inverse(),
                "micron"
            )
            set_transformation(
                sdata_main[f"boundaries_{seg_method}"],
                Identity(),
                "pixel"
            )
    else:
        set_transformation(
            sdata_main[f"boundaries_{seg_method}"],
            Identity(),
            "micron"
        )
        set_transformation(
            sdata_main[f"boundaries_{seg_method}"],
            transformation_to_pixel,
            "pixel"
        )
    return


def transform_adata(
    sdata_main: sd.SpatialData, seg_method: str, data_path: str
) -> None:
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

def get_2D_boundaries(method: str, org_sdata:sd.SpatialData, sdata:sd.SpatialData, transformation:np.ndarray, boundary_key:str) -> None:
    if method.startswith("vpt_3D"):
        try:
            bound = org_sdata[boundary_key][["cell_id", "geometry"]].dissolve(by="cell_id")
        except ValueError:
            new = org_sdata[boundary_key][["cell_id", "geometry"]]
            new.index.rename(None, inplace=True)
            bound = new.dissolve(by="cell_id")
            del new
        sdata[f"boundaries_{method}"] = ShapesModel.parse(bound)
        set_transformation(
            sdata[f"boundaries_{method}"],
            Affine(transformation, input_axes=("x", "y"), output_axes=("x", "y")),
            to_coordinate_system="global",
        )
    else:
        sdata[f"boundaries_{method}"] = ShapesModel.parse(org_sdata[boundary_key])
    if any([method.startswith(x) for x in image_based]):
        if method == "Cellpose_1_Merlin":
            set_transformation(
                sdata[f"boundaries_{method}"],
                Identity(),
                to_coordinate_system="micron",
            )
        else:
            set_transformation(
                sdata[f"boundaries_{method}"],
                Affine(transformation, input_axes=("x", "y"), output_axes=("x", "y")).inverse(),
                to_coordinate_system="micron",
            )
    else:
        set_transformation(
            sdata[f"boundaries_{method}"],
            Identity(),
            to_coordinate_system="micron",
        )


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
    top_n_factors: int = 3,
    n_ficture: int = 21,
    logger: logging.Logger = None,
    factors: Optional[List[int]] = None,
) -> Dict[str, Union[np.ndarray, List[int]]]:
    """Generate ficture images stack and other ficture information.

    Args:
        data_path: Path to merscope data
        results_path: path to Segmentation folder
        top_n_factors: only consider top n factors for ficture picture
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
    # assume, that transform[0,1], transform[1,0] = 0
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
        unique_factors = unique_factors.union(set(np.unique(ficture_pixels[f"K{i}"])))
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
            image_stack = create_factor_level_image(
                ficture_pixels, factor, DAPI_shape, top_n_factors
            )
        else:
            image_stack = np.concatenate(
                (
                    image_stack,
                    create_factor_level_image(
                        ficture_pixels, factor, DAPI_shape, top_n_factors
                    ),
                ),
                axis=0,
                dtype=np.uint16,
            )
    return {"images": image_stack, "factors": unique_factors}


def compute_cell_morphology(
    sdata, add_to_adata=False, z_spacing=1.5, verbose=True, return_results=False
):
    """Compute 2D or 3D morphology metrics from shape boundaries in a SpatialData object.

    Metrics per cell include:
    - dimensionality: '2D' or '3D', depending on whether the cell spans multiple Z planes.
    - area: Total surface area (2D: shape area; 3D: approximated surface area).
    - volume_sum: Naive volume by summing areas and multiplying with z_spacing.
    - volume_trapz: Trapezoidal volume integration across z-planes (3D only).
    - volume_final: Final volume estimate (trapz if 3D, volume_sum if 2D).
    - num_z_planes: Number of Z planes a cell spans.
    - size_normalized: Edge length of a square (2D) or cube (3D) with same area or volume as the cell.
    - surface_to_volume_ratio: Ratio of surface area to volume.
    - sphericity: Shape compactness; 1 for a perfect sphere (3D only).
    - solidity: Volume compared to its convex hull volume (compactness).
    - elongation: PCA-based anisotropy estimate (0 = isotropic, 1 = elongated).
    """
    results_dict = {}

    for boundaries_name, boundaries in sdata.shapes.items():
        if not boundaries_name.startswith("boundaries_"):
            continue

        dataset_name = boundaries_name.replace("boundaries_", "")
        table_name = f"adata_{dataset_name}"

        if table_name not in sdata.tables:
            warnings.warn(f"Missing adata {table_name}. Skipping {boundaries_name}.")
            continue

        if boundaries.empty or "geometry" not in boundaries:
            warnings.warn(f"Empty/invalid geometry in {boundaries_name}. Skipping.")
            continue

        # Parse geometry from strings if needed
        boundaries = boundaries.copy()
        boundaries["geometry"] = boundaries["geometry"].apply(
            lambda x: Polygon(ast.literal_eval(x)) if isinstance(x, str) else x
        )

        if "EntityID" not in boundaries:
            warnings.warn(f"Missing EntityID in {boundaries_name}. Skipping.")
            continue

        global_z_min = boundaries["ZIndex"].min()
        global_z_max = boundaries["ZIndex"].max()

        grouped = boundaries.groupby("EntityID")
        morphology_rows = []

        for entity_id, group in grouped:
            try:
                is_entity_3d = "ZIndex" in group and group["ZIndex"].nunique() > 1
                polygons = group["geometry"].tolist()

                if not polygons:
                    continue

                if is_entity_3d:
                    morphology_data = _compute_3d_metrics(
                        group,
                        polygons,
                        z_spacing,
                        global_z_min=global_z_min,
                        global_z_max=global_z_max,
                        verbose=verbose,
                    )
                else:
                    morphology_data = _compute_2d_metrics(polygons[0], z_spacing)

                morphology_data["EntityID"] = entity_id
                morphology_rows.append(morphology_data)
            except Exception as e:
                warnings.warn(f"Failed to process entity {entity_id}: {str(e)}")
                continue

        if morphology_rows:
            df_morph = pd.DataFrame(morphology_rows).set_index("EntityID")

            if add_to_adata:
                adata = sdata.tables[table_name]
                adata_ids = adata.obs.index.astype(str)
                common = adata_ids.intersection(df_morph.index.astype(str))
                for col in df_morph.columns:
                    adata.obs[f"cell_{col}"] = np.nan
                    adata.obs.loc[common, f"cell_{col}"] = df_morph.loc[common, col]

            if verbose:
                dim_counts = df_morph["dimensionality"].value_counts().to_dict()
                print(f"Added morphology metrics for {dataset_name}: {dim_counts}")

            if return_results:
                results_dict[dataset_name] = df_morph

    return results_dict if return_results else None


def _compute_2d_metrics(geom, z_spacing: float):
    """Compute 2D morphology metrics with improved robustness."""
    if geom.geom_type == "MultiPolygon":
        geom = unary_union(geom)
    if geom.geom_type == "MultiPolygon":
        geom = max(geom.geoms, key=lambda p: p.area)

    if geom.is_empty or geom.geom_type != "Polygon":
        return {}

    area = geom.area
    perimeter = geom.length

    metrics = {
        "dimensionality": "2D",
        "area": area,
        "volume_sum": area * z_spacing,
        "volume_final": area * z_spacing,
        "num_z_planes": 1,
        "size_normalized": np.sqrt(area),
        "surface_to_volume_ratio": perimeter / area if area > 0 else np.nan,
    }

    metrics["sphericity"] = (
        4 * PI * area / (perimeter**2) if perimeter > 0 else np.nan
    )

    if not geom.is_valid:
        geom = geom.buffer(0) #assume no topology error occurs
    convex_hull = geom.convex_hull
    metrics["solidity"] = (
        area / convex_hull.area if convex_hull.area > 0 else np.nan
    )

    hull_points = np.array(convex_hull.exterior.coords[:-1])
    if len(hull_points) >= 3:
        hull_points -= hull_points.mean(axis=0)
        hull_points /= hull_points.std(axis=0, ddof=0)
        cov = (hull_points.T @ hull_points) / hull_points.shape[0]
        eigenvalues  = np.linalg.eigvalsh(cov)[::-1]
        metrics["elongation"] = (
            1 - np.sqrt(eigenvalues[1] / eigenvalues[0])
            if eigenvalues[0] > 0
            else np.nan
        )
    else:
        metrics["elongation"] = np.nan

    return metrics


def _compute_trapz_with_conditional_caps(areas, z_indices, z_spacing, z_min, z_max):
    """Compute trapezoidal cell volume + add half caps only if top/bottom areas are at imaging boundaries (hence, assuming cell is truncated)."""
    if len(areas) == 1:
        return areas[0] * z_spacing
    volume = np.trapezoid(areas, dx=z_spacing)
    if z_indices[0] == z_min:
        volume += 0.5 * z_spacing * areas[0]
    if z_indices[-1] == z_max:
        volume += 0.5 * z_spacing * areas[-1]
    return volume


def _compute_3d_metrics(
    group, z_spacing, global_z_min, global_z_max, ZIndex: str="ZIndex"
):
    """Compute 3D morphology metrics with robust and dual volume estimation."""
    try:
        group = group.sort_values(ZIndex)
        z_indices = group[ZIndex].to_numpy()
        z_um = z_indices * z_spacing
        polygons = group["geometry"].tolist()

        if z_indices.shape[0] < 1:
            return _compute_2d_metrics(polygons, z_spacing)

        assert len(polygons) == len(z_indices)
        assert np.all(np.diff(z_indices) >= 0), "ZIndex must be non-decreasing"

        areas = np.fromiter(
            (
                p.area for p in polygons #Multipolygons bereits gehandled
            ), dtype=float
        )

        volume_sum = np.sum(areas) * z_spacing
        volume_trapz = _compute_trapz_with_conditional_caps(
            areas, z_indices, z_spacing, global_z_min, global_z_max
        )

        perimeters = np.fromiter(
            (
                p.length for p in polygons #Multipolygons bereits gehandled
            ), dtype=float
        )

        total_height = (len(areas) - 1) * z_spacing
        lateral_surface = (
            np.mean(perimeters) * total_height if len(perimeters) > 1 else 0
        )
        top_bottom_surface = areas[0] + areas[-1] if len(areas) > 0 else 0
        surface_area = lateral_surface + top_bottom_surface

        metrics = {
            "dimensionality": "3D",
            "area": surface_area,
            "volume_sum": volume_sum,
            "volume_trapz": volume_trapz,
            "volume_final": volume_trapz,
            "num_z_planes": len(areas),
            "size_normalized": np.cbrt(volume_trapz) if volume_trapz > 0 else np.nan,
            "surface_to_volume_ratio": surface_area / volume_trapz
            if volume_trapz > 0
            else np.nan,
            "sphericity": (
                (PI ** (1 / 3)) * (6 * volume_trapz) ** (2 / 3) / surface_area
                if volume_trapz > 0 and surface_area > 0
                else np.nan
            ),
        }

        thin_stack = len(polygons) < 2 or (np.ptp(z_um) < 3 * z_spacing)

        all_points = []
        for i, polygon in enumerate(polygons):
            if polygon.is_empty or not polygon.is_valid:
                continue
            if polygon.geom_type == "Polygon":
                coords = np.array(polygon.exterior.coords[:-1])
            else:
                coords = np.vstack([np.asarray(part.exterior.coords)[:-1] for part in polygon.geoms]) \
                                        if len(polygon.geoms) else np.empty((0, 2))
            if coords.shape[0] < 3:
                continue
            z_coords = np.full((coords.shape[0], 1), z_um[i])
            all_points.append(np.hstack([coords, z_coords]))

        if not all_points:
            metrics.update({"solidity": np.nan, "elongation": np.nan})
            return metrics

        all_points = np.vstack(all_points)
        all_points -= all_points.mean(axis=0, keepdims=True)

        solidity = np.nan
        if not thin_stack:
            if all_points.shape[0] >= 4:
                try:
                    hull = ConvexHull(all_points, qhull_options="QJ")  # jitter coplanar cases
                    hv = hull.volume
                    solidity = (volume_trapz / hv) if hv > 0 else np.nan
                except:
                    # fall back to deduped points (expensive; only if needed)
                    up = np.unique(all_points, axis=0)
                    if up.shape[0] >= 4:
                        try:
                            hull = ConvexHull(up, qhull_options="QJ")
                            hv = hull.volume
                            solidity = (volume_trapz / hv) if hv > 0 else np.nan
                        except:
                            solidity = np.nan
        metrics["solidity"] = solidity

        cov = (all_points.T @ all_points) / all_points.shape[0]
        eigenvalues = np.linalg.eigvalsh(cov)[::-1]

        if eigenvalues[0] <= 0:
            elongation = np.nan
        elif eigenvalues[2] > 1e-6:
            elongation = 1 - np.sqrt(eigenvalues[2] / eigenvalues[0])
        elif eigenvalues[1] > 1e-6:
            elongation = 1 - np.sqrt(eigenvalues[1] / eigenvalues[0])
        else:
            elongation = 1.0
        metrics["elongation"] = elongation

        return metrics

    except Exception as e:
        warnings.warn(f"Failed to compute 3D metrics: {str(e)}")
        return {}
