import sys
import os
import logging
from os.path import join
from pathlib import Path
from collections import defaultdict
from re import split

from anndata import AnnData
from spatialdata import read_zarr

sys.path.insert(1, join(str(Path(__file__).parent.parent.resolve()), "src"))
from ficture_utils import create_factor_level_image, parse_metadata
from metrics.ficture_intensities import ficture_intensities
from sdata_utils import build_shapes, add_cell_type_annotation, add_ficture, transform_adata, prepare_ficture, update_element

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def check_cell_type_availability(adata: AnnData, method_path: str)->bool:
    if "cell_type_annotation" in os.listdir(method_path):
        if "cell_type_final" in list(adata.obs.columns):
            return False
        else:
            return True
    else:
        return False

def check_ficture_availability(adata: AnnData, sdata_path: str, n_ficture:int, var=False)->bool:
    if "Ficture" not in os.listdir(join(sdata_path, "results")):
        return False

    if not var:
        if f"ficture{n_ficture}_mean" in list(adata.obsm_keys()):
            return False
    else:
        if set([f"ficture{n_ficture}_mean", f"ficture{n_ficture}_var"]) <= set(adata.obsm_keys()) :
            return False

    ficture_path = join(sdata_path, "results", "Ficture", "output")
    for file in os.listdir(ficture_path):
        if n_ficture == int(split(r"\.|F", file)[1]):
            return True
    return False

sample = sys.argv[1]  # sample name
data_path = sys.argv[2]  # merscope data
zmode = sys.argv[3]
data_dir = sys.argv[4]  # base directory
if len(sys.argv) > 5:
    n_ficture = int(sys.argv[5])
    if len(sys.argv) > 6:
        var = sys.argv[6]
    else:
        var = False
else:
    n_ficture = 21
    var = False

sdata_path = join(data_dir, "samples", sample)
sdata_main = read_zarr(join(sdata_path, "sdata_z3.zarr"))

logging.info("Starting task collection")
seg_methods = [
    method
    for method in os.listdir(join(sdata_path, "results"))
    if os.path.isdir(join(sdata_path, "results", method, "sdata.zarr"))
]

print("Compile task list")
tasks_collection = defaultdict(list)

ficture_arguments = prepare_ficture(data_path, sdata_path, n_ficture)
if not ficture_arguments:
    logging.warning(f"No ficture output found.")

for method in seg_methods:
    if f"boundaries_{method}" not in sdata_main.shapes.keys():
        tasks_collection[f"boundaries_{method}"].append(f"shapes")
    if f"adata_{method}" not in sdata_main.tables.keys():
        tasks_collection[f"adata_{method}"].append('adata')
    else:
        adata = sdata_main.tables[f"adata_{method}"]
        method_path = os.path.join(sdata_path, method)

        if not os.path.exists(join(method_path, "cell_type_annotation")):
            logging.warning(f"No cell type annotation found on disk for method '{method}'")
        elif "cell_type_final" not in adata.obs.columns:
            tasks_collection[f"adata_{method}"].append("cell_types")

        if check_ficture_availability(adata, sdata_path, n_ficture, var=var):
            tasks_collection[f"adata_{method}"].append("ficture")

for key, upd in tasks_collection.items():
    logging.info(f"{key} requires updates: {upd}")

for method, tasks in tasks_collection.items():
    logging.info(f"Starting updates for '{method}'")
    sdata = read_zarr(join(sdata_path, "results", method, "sdata.zarr"))
    for task in tasks:
        logging.info(f"Applying task '{task}' for '{method}'")
        if task == "shapes":
            sdata_main = build_shapes(sdata, sdata_main, method, write_to_disk=False)
        if task == "cell_types":
            sdata_main = add_cell_type_annotation(sdata_main, sdata_path, method, write_to_disk=False)
        elif task == "ficture":
            sdata_main = add_ficture(sdata, sdata_main, method, ..., n_ficture, var, write_to_disk=False)
        elif task == "adata":
            sdata_main[f"adata_{method}"] = sdata[list(sdata.tables.keys())[0]].copy()
            transform_adata(sdata_main, method, data_path=data_path)

            if "cell_type_annotation" in os.listdir(
                    join(sdata_path, "results", method)
            ):  # TODO: add automatic cell type annotation
                sdata_main = add_cell_type_annotation(sdata_main, sdata_path, method, write_to_disk=False)
            else:
                logging.warning(f"No cell type annotation found for '{method}'. Skipping annotation.")

            if len(ficture_arguments) > 0:
                sdata_main = add_ficture(sdata, sdata_main, method, ficture_arguments, n_ficture, var,
                                         write_to_disk=False)
    if "shapes" in tasks:
        update_element(sdata_main, f"boundaries_{method}")
        logging.info(f"Completed shape update for '{method}'")
    elif len(set(["cell_type_annotation", "ficture", "adata"]) & set(tasks)):
        update_element(sdata_main, f"adata_{method}")
        logging.info(f"Completed adata update for '{method}'")
