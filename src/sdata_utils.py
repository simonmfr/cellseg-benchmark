import os
import glob
import spatialdata as sd
import spatialdata_io

def process_merscope(sample_name, data_dir, sample_paths, zmode):
    """Load and save a MERSCOPE sample as sdata with specified z_layers configuration."""
    if zmode not in {"z3", "3d"}:
        raise ValueError(f"Invalid zmode: {zmode}")
    sdata_file = os.path.join(data_dir, "samples", sample_name, f"sdata_{zmode}.zarr")
    if os.path.exists(sdata_file):
        print(f"Skipping {sample_name}: {zmode} file already exists")
        return   
    sdata = spatialdata_io.merscope(
        sample_paths[sample_name],
        z_layers = 3 if zmode == "z3" else range(7),
        backend=None,
        cells_boundaries=False,
        cells_table=False,
        mosaic_images=True,
        transcripts=True,
        slide_name="_".join(sample_name.split("_")[:2]),
        region_name=sample_name.split("_")[2],
    )
    os.makedirs(os.path.dirname(sdata_file), exist_ok=True)
    sdata.write(sdata_file, overwrite=False)


def process_merlin_segmentation(sample_name, sample_paths, sdata_main, write_to_disk=True):
    """
    Process Merlin-specific segmentation data and add it to the main spatial data object.
    
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
    if f"boundaries_{seg_method}" not in sdata_main or f"adata_{seg_method}" not in sdata_main:
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
         
    return sdata_main


def integrate_segmentation_data(data_dir, sample_name, seg_methods, sdata_main, write_to_disk=True):
    """
    Integrate segmentation data from multiple methods into the main spatial data object.
    
    Args:
        data_dir: Base directory for data
        sample_name: Name of the sample
        seg_methods: List of segmentation methods to process
        sdata_main: Main spatial data object to update
        write_to_disk: Whether to write elements to disk immediately
    
    Returns:
        Updated sdata_main object
    """    
    for seg_method in seg_methods:
        seg_path = os.path.join(data_dir, "samples", sample_name, 
                               "results", seg_method, "sdata.zarr")
        if not os.path.exists(seg_path):
            continue
            
        if f"boundaries_{seg_method}" not in sdata_main or f"adata_{seg_method}" not in sdata_main:
            sdata = sd.read_zarr(seg_path)
        
        # Handle boundaries
        if f"boundaries_{seg_method}" not in sdata_main:
            base_dir = os.path.join(seg_path, "shapes")
            boundary_files = [os.path.basename(f) for f in glob.glob(os.path.join(base_dir, "*_boundaries"))]
            
            if len(boundary_files) == 1:
                sdata_main[f"boundaries_{seg_method}"] = sdata[boundary_files[0]]
                if write_to_disk:
                    sdata_main.write_element(f"boundaries_{seg_method}")
            elif len(boundary_files) > 1:
                print(f"Multiple *boundaries files found for {seg_method}. Skipping.")
        else:
            print(f"Skipping {seg_method} as boundaries_{seg_method} exist already.")
        
        # Handle tables
        if f"adata_{seg_method}" not in sdata_main:
            base_dir = os.path.join(seg_path, "tables")
            table_files = [os.path.basename(f) for f in glob.glob(os.path.join(base_dir, "table"))]
            
            if len(table_files) == 1:
                sdata_main[f"adata_{seg_method}"] = sdata[table_files[0]]
                if write_to_disk:
                    sdata_main.write_element(f"adata_{seg_method}")
            elif len(table_files) > 1:
                print(f"Multiple table files found for {seg_method}. Skipping.")
            else:
                print(f"Table file missing for {seg_method}. Skipping.")
        else:
            print(f"Skipping {seg_method} as adata_{seg_method} exist already.")
    
    return sdata_main