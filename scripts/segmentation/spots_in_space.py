import numpy as np
import os
from tqdm import tqdm

import argparse
from pathlib import Path
import dask
dask.config.set({"dataframe.query-planning": False})
import sis

parser = argparse.ArgumentParser(description="Run sis cellpose segmentation on a MERSCOPE slide.")
parser.add_argument("data_path", type=Path, help="Path to merfish output folder.")
parser.add_argument("save_path", type=Path, help="Path to sis output folder.")
parser.add_argument("model_path", type=Path, help="Path to cellpose model.")
parser.add_argument("staining", help="Cytoplasm channel (e.g. total_mrna, PolyT).")
args = parser.parse_args()

def main(data_path, save_path, model_path, staining):
    """Spots-in-space segmentation using cellpose 3D, writes polygons and cell-by-gene."""
    save_path.mkdir(exist_ok=True, parents=True)
    cache_file = save_path / "detected_transcripts.csv.npz"
    st = sis.SpotTable.load_merscope(
        csv_file=data_path / "detected_transcripts.csv",
        cache_file=cache_file,
        image_path=data_path / "images",
    )
    
    # test on subset
    #subrgn = ((12000, 12200), (2500, 2700))
    #st = st.get_subregion(xlim=subrgn[0], ylim=subrgn[1])
    
    seg_opts = {
        "cellpose_model": str(model_path),
        "cellpose_gpu": "auto",
        "px_size": 0.108,
        "cell_dia": 10,
        "z_plane_thickness": 1.5,
        "images": {
            "cyto":   {"channel": staining, "n_planes": 7},
            "nuclei": {"channel": "DAPI"},
        },
        "cellpose_options": {"batch_size": 8, "min_size": 5000},
    }
    
    seg_custom_3d = sis.segmentation.CellposeSegmentationMethod(seg_opts)
    tiles = st.grid_tiles(max_tile_size=350, overlap=30, min_transcripts=1000) # default used by Allen Institute
    for i, tile in enumerate(tqdm(tiles)):
        result = seg_custom_3d.run(tile)
        np.save(save_path / f'cell_ids_{i}.npy', result.cell_ids)
    
    truncated_meta = {
                'seg_method': f"sis_dapi_{staining}",
                'seg_opts': seg_opts,
                #'polygon_opts': polygon_opts
                }
    
    seg_st = sis.spot_table.SegmentedSpotTable(
                    spot_table=st, 
                    cell_ids=np.empty(len(st), dtype=int),
                    seg_metadata=truncated_meta,
                    )
    
    merge_results = []
    skipped = []
    for i, tile in enumerate(tiles):
        cell_id_file = save_path / f'cell_ids_{i}.npy'
        if not os.path.exists(cell_id_file):
            print(f"Skipping tile {i} : no cell ID file generated")
            skipped.append(i)
            continue
    
        tile_cids = np.load(cell_id_file)
        tile = sis.spot_table.SegmentedSpotTable(tile, tile_cids)
    
        tiles[i] = tile
    
    _ = seg_st.set_cell_ids_from_tiles(tiles, padding=seg_opts['cell_dia'] / 2)
    
    seg_st.calculate_cell_polygons()
    seg_st.save_cell_polygons(save_path / "cell_polygons.geojson")
    
    seg_st.generate_cell_labels(prefix=None, suffix=None)
    cell_by_gene = seg_st.cell_by_gene_anndata(x_format='sparse')
    seg_st.save_anndata(save_path / "cell_by_gene.h5ad", cell_by_gene)
    
    cache_file.unlink(missing_ok=True)
    for f in save_path.glob("cell_ids_*.npy"):
        f.unlink()
    print("Done.")
    
if __name__ == "__main__":
    main(args.data_path, args.save_path, args.model_path, args.staining)