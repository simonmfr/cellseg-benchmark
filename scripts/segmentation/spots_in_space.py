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
    try:
        st = sis.SpotTable.load_merscope(
            csv_file=data_path / "detected_transcripts.csv",
            cache_file=cache_file,
            image_path=data_path / "images",
        )
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
        seg = sis.segmentation.CellposeSegmentationMethod(seg_opts)
        result = seg.run(st)
        seg_st = result.spot_table()
        seg_st.calculate_cell_polygons()
        seg_st.save_cell_polygons(save_path / "cell_polygons.geojson")
        seg_st.generate_cell_labels()
        seg_st.save_anndata(save_path / "cell_by_gene.h5ad", x_format="sparse")
        print("done")
    finally:
        cache_file.unlink(missing_ok=True)


if __name__ == "__main__":
    main(args.data_path, args.save_path, args.model_path, args.staining)