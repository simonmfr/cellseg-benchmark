# User guide

## api_baysor.py
Arguments:
1) confidence
2) base segmentation: full name, i.e. Cellpose_1_DAPI_PolyT
3) sample name

Usage:
```
python api_baysor.py confidence base_segmentation sample_name
```

## api_cellpose_1.py
Arguments:
1) data path
2) save path
3) staining: cytoplasma staining for segmentation

Usage:
```
python api_cellpose_1.py data_path save_path staining
```

## api_cellpose_2.py
Arguments:
1) data path
2) save path
3) staining: cytoplasma staining for segmentation

Usage:
```
python api_cellpose_2.py data_path save_path staining
```

## api_nuclei.py
Arguments:
1) data path
2) save path

Usage:
```
python api_nuclei.py data_path save_path
```

## rastered_segmentation.py
Arguments:
1) data path
2) save path
3) width: width of raster
4) overlap: Overlap of cells
5) intensity ratio: parameter for filter of cells

Usage:
```
python rastered_segmentation.py data_path save_path width overlap intensity_ratio
```

## voronoi_segmentation.py
Arguments:
1) data path
2) save path
3) max number of cells

Usage:
```
python voronoi_segmentation.py data_path save_path max_cells
```