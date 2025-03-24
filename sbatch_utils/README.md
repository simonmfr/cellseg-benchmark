# User Guide
Creating sbatch scripts for different methods

## [annotation_script_creation.py](annotation_script_creation.py)
Arguments:
1) `sample`: name of sample

Usage:
```
python annotation_script_creation.py sample
```

## [baysor_script_creation.py](baysor_script_creation.py)
Arguments:
1) `staining`: Staining of the prior cellpose segmentation
2) `cellpose_version`: Cellpose version for prior segmentation
3) `confidence`: confidence in prior segmentation

Usage:
```
python baysor_script_creation.py staining cellpose_version confidence
```

## [cellpose1_script_creation.py](cellpose1_script_creation.py)
Arguments:
1) `sample`: name of sample

Usage:
```
python cellpose1_script_creation.py sample
```

## [cellpose2_script_creation.py](cellpose2_script_creation.py)
Arguments:
1) `sample`: name of sample

Usage:
```
python cellpose2_script_creation.py sample
```

## [rastered_scripts.py](rastered_scripts.py)
Arguments:
1) `length`: length of the sides of the squares
2) `unit`: either in microns or pixels
3) `overlap`: Overlap between squares
4) `intensity_ratio`: between 1 and 0. Otherwise defaults to 0.1

Usage:
```
python rastered_scripts.py length unit overlap intensity_ratio
```