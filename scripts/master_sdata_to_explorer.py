from operator import contains
from os.path import join
from sys import argv

from anndata import concat as adconcat
from pandas import Categorical
from pandas import concat as pdconcat
from sopa.io.explorer import write
from spatialdata import read_zarr
from spatialdata.models import ShapesModel


sdata = read_zarr(join(argv[1], "sdata_z3.zarr"))
methods = argv[2:]

assert all([method in sdata.tables.keys() for method in methods]), "Not all keys are contained in master sdata."

# cell_id needs to be in the adata for processing reasons
# keep the name of the method
adatas_name = [
    "_".join(key.split("_")[1:])
    for key in list(sdata.tables.keys())
    if "cell_id" in sdata[key].obs.keys() and key in methods
]

adatas_shapes = ["boundaries_" + shape for shape in adatas_name]
shapes_collected = pdconcat(
    [sdata[key] for key in list(sdata.shapes.keys()) if key not in adatas_shapes]
)

adatas_tables = ["adata_" + table for table in adatas_name]
adatas_tables = []
for key in list(sdata.tables.keys()):
    if key in adatas_tables:
        adata = sdata[key]
        adata.obs["method"] = "_".join(key.split("_")[1:])  # name method explicitly
        adata.obs["method"] = Categorical(adata.obs["method"])
        adatas_tables.append(adata)
tables_collected = adconcat(adatas_tables, uns_merge="first")

tables_collected.obs["region_new"] = "shapes_collected"
tables_collected.uns["spatialdata_attrs"]["region"] = "shapes_collected"
tables_collected.uns["spatialdata_attrs"]["region_key"] = "region_new"

sdata["shapes_collected"] = ShapesModel.parse(shapes_collected)
sdata["tables_collected"] = tables_collected

write(
    join(argv[1], "sdata_z3.explorer"),
    sdata,
    table_key="tables_collected",
    shapes_key="shapes_collected",
    gene_column="gene",
    save_h5ad=True,
)
