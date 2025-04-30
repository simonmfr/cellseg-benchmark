# User guide

## [api_baysor.py](api_baysor.py)
Arguments:
1) `base segmentation`: full name, i.e. Cellpose_1_DAPI_PolyT
2) `confidence`: confidence of prior segmentation
3) `sample_name`: name of sample

Usage:
```
python api_baysor.py base_segmentation confidence sample_name
```

## [api_cellpose_1.py](api_cellpose_1.py)
Arguments:
1) `data_path`: Path to Merscope output
2) `save_path`: Path to results folder
3) `staining`: cytoplasma staining for segmentation

Usage:
```
python api_cellpose_1.py data_path save_path staining
```

## [api_cellpose_2.py](api_cellpose_2.py)
Arguments:
1) `data_path`: Path to Merscope output
2) `save_path`: Path to results folder
3) `staining`: cytoplasma staining for segmentation

Usage:
```
python api_cellpose_2.py data_path save_path staining
```

## [api_nuclei.py](api_nuclei.py)
Arguments:
1) `data_path`: Path to Merscope output
2) `save_path`: Path to results folder

Usage:
```
python api_nuclei.py data_path save_path
```

## [rastered_segmentation.py](rastered_segmentation.py)
Arguments:
1) `data_path`: Path to Merscope output
2) `save_path`: Path to results folder
3) `width`: width of raster
4) `overlap`: Overlap of cells
5) `intensity_ratio`: parameter for filter of cells

Usage:
```
python rastered_segmentation.py data_path save_path width overlap intensity_ratio
```

## [voronoi_segmentation.py](voronoi_segmentation.py)
Arguments:
1) `data_path`: Path to Merscope output
2) `save_path`: Path to results folder
3) `max_cells`: Max number of cells for segmentation

Usage:
```
python voronoi_segmentation.py data_path save_path max_cells
```

## [run_mapmycells.py](run_mapmycells.py)
Run automated cell type annotation using MapMyCells (MMC) from Allen Institute and revise output.

1) Run MMC on mouse brain reference atlas (ABCAtlas)
2) QCs MMC output, defined low-quality mappings as "Undefined", and mixed cell type identities as "Mixed", then merged labels into adata
3) To de-noise data, run Leiden clustering and assign cell types to clusters
4) Revise annotation using marker gene scores from reference atlas
5) UMAP and Spatial plots
6) Export results as csv
   
Arguments:
1. `sample_name`: Name of the sample to process (eg foxf2_s2_r1)
2. `method_name`: Name of the segmentation method (eg Proseg)
3. `data_directory`: Root of analysis folder (cellseg-benchmark)
5. `mad_factor`: Factor for median absolute deviation filtering (defaults to 3 if not provided or if value < 1)

Usage:
```bash
python run_mapmycells.py sample_name method_name data_directory [mad_factor]
```

## [transcript_tif.py](transcript_tif.py)
Arguments:
1) `data_path`: Path to Merscope output

Usage:
```
python transcript_tif.py data_path
```

## [master_sdata_to_explorer.py](master_sdata_to_explorer.py)
Convert multiple shapes and boundaries into one explorer.
Arguments:
1) `data_path`: Path to sdata_z3.zarr
2) `method_1`: Method to add to explorer
3) `method_2`: Method to add to explorer
4) `...`: More methods to add to explorer

Usage:
```
python master_sdata_to_explorer.py data_path method_1 method_2 ...
```