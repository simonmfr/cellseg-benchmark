from spatialdata import read_zarr
from sopa.io.explorer import write
from spatialdata.models import ShapesModel
from pandas import Categorical
from pandas import concat as pdconcat
from anndata import concat as adconcat
from sys import argv
from os.path import join

sdata = read_zarr(join(argv[1], "sdata_z3.zarr"))

# cell_id needs to be in the adata for processing reasons
# keep the name of the method
bad_adatas = ["_".join(key.split("_")[1:]) for key in list(sdata.tables.keys()) if "cell_id" not in sdata[key].obs.keys()]

bad_adatas_shapes = ["boundaries_"+shape for shape in bad_adatas]
shapes_collected = pdconcat([sdata[key] for key in list(sdata.shapes.keys()) if key not in bad_adatas_shapes])

bad_adatas_tables = ["adata_"+table for table in bad_adatas]
adatas_tables = []
for key in list(sdata.tables.keys()):
    if key not in bad_adatas_tables:
        adata = sdata[key]
        adata.obs['method'] = "_".join(key.split("_")[1:]) #name method explicitly
        adata.obs['method'] = Categorical(adata.obs['method'])
        adatas_tables.append(adata)
tables_collected = adconcat(adatas_tables, uns_merge="first")

tables_collected.obs["region_new"] = "shapes_collected"
tables_collected.uns["spatialdata_attrs"]['region'] = "shapes_collected"
tables_collected.uns["spatialdata_attrs"]['region_key'] = "region_new"

sdata["shapes_collected"] = ShapesModel.parse(shapes_collected)
sdata["tables_collected"] = tables_collected

write(join(argv[1], "sdata_z3.explorer"), sdata, table_key='tables_collected', shapes_key='shapes_collected', gene_column="gene", save_h5ad=True)
