from setuptools import find_packages, setup

setup(
    name="cellseg_benchmark",
    version="0.1.0",
    description="Tools for assessing the quality of spatial segmentations with focus on vascular cells in the brain.",
    author="IDS Munich",
    packages=find_packages(exclude=["notebooks", "scripts"]),
    install_requires=[
        "sopa==2.0.6",
        "spatialdata",
        "spatialdata-io>=0.1.5",
        "scanpy",
        "scikit-learn",
        "scipy",
        "seaborn",
        "numpy",
        "pandas",
        "shapely",
        "dask",
        "geopandas",
        "xarray",
        "tifffile",
        "tqdm",
        "anndata",
        "scikit-image",
        "plotly",
        "ovrlpy"
    ],
    python_requires=">=3.10",
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
)
