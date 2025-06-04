import logging
import os
import sys
import warnings
from collections import defaultdict
from os.path import exists, join
from re import split

from anndata import AnnData
from spatialdata import read_zarr

from cellseg_benchmark.sdata_utils import (
    add_cell_type_annotation,
    add_ficture,
    build_shapes,
    calculate_volume,
    prepare_ficture,
    transform_adata,
    update_element,
)

logger = logging.getLogger("shape_mapping")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

warnings.filterwarnings("ignore")


def check_ficture_availability(
    adata: AnnData, sdata_path: str, n_ficture: int, var=False
) -> bool:
    """Check if ficture is available."""
    if "Ficture" not in os.listdir(join(sdata_path, "results")):
        return False

    if not var:
        if f"ficture{n_ficture}_mean" in list(adata.obsm_keys()):
            return False
    else:
        if set([f"ficture{n_ficture}_mean", f"ficture{n_ficture}_variance"]) <= set(
            adata.obsm_keys()
        ):
            return False

    ficture_path = join(sdata_path, "results", "Ficture", "output")
    for file in os.listdir(ficture_path):
        if n_ficture == int(split(r"\.|F", file)[1]):
            return True
    return False


def clean_tasks(tasks):
    """Clean the tasks."""
    tasks_tmp = set(tasks)
    if "adata" in tasks:
        tasks_tmp.discard("cell_types")
        tasks_tmp.discard("volume")
        tasks_tmp.discard("ficture")
    return tasks_tmp


sample = sys.argv[1]  # sample name
data_path = sys.argv[2]  # merscope data
zmode = sys.argv[3]
data_dir = sys.argv[4]  # base directory
if len(sys.argv) > 5:
    n_ficture = int(sys.argv[5])
else:
    n_ficture = 21
if len(sys.argv) > 6:
    var = True
else:
    var = False
if len(sys.argv) > 7:
    if sys.argv[7] == "--force":
        forcing = sys.argv[8:]
    else:
        forcing = []
else:
    forcing = []

sdata_path = join(data_dir, "samples", sample)
sdata_main = read_zarr(join(sdata_path, "sdata_z3.zarr"))

logger.info("Starting task collection")
seg_methods = [
    method
    for method in os.listdir(join(sdata_path, "results"))
    if os.path.isdir(join(sdata_path, "results", method, "sdata.zarr"))
]

tasks_collection = defaultdict(list)

if "Ficture" not in os.listdir(join(sdata_path, "results")):
    logger.warning("No ficture output found.")

ficture_flag = False
for method in seg_methods:
    if forcing != []:
        tasks_collection[method].extend(forcing)
    if f"boundaries_{method}" not in sdata_main.shapes.keys():
        tasks_collection[method].append("shapes")
    if f"adata_{method}" not in sdata_main.tables.keys():
        if exists(join(sdata_path, "results", method, "sdata.zarr", "tables")):
            tasks_collection[method].append("adata")
    else:
        adata = sdata_main.tables[f"adata_{method}"]
        method_path = os.path.join(sdata_path, "results", method)

        if not os.path.exists(
            join(method_path, "cell_type_annotation", "adata_obs_annotated.csv")
        ):
            logger.warning(
                f"No cell type annotation found on disk for method '{method}'"
            )
        else:
            if not any(["cell_type" in x for x in adata.obs.columns]):
                tasks_collection[method].append("cell_types")

        if "volume" not in adata.obs.columns:
            tasks_collection[method].append("volume")

        if check_ficture_availability(adata, sdata_path, n_ficture, var=var):
            tasks_collection[method].append("ficture")
            ficture_flag = True

for key, upd in tasks_collection.items():
    logger.info(f"{key} requires updates: {upd}")

if ficture_flag:
    ficture_arguments = prepare_ficture(data_path, sdata_path, n_ficture)
else:
    ficture_arguments = None

for method, tasks in tasks_collection.items():
    logger.info(f"Starting updates for '{method}'")
    tasks = clean_tasks(tasks)
    sdata = read_zarr(join(sdata_path, "results", method, "sdata.zarr"))
    for task in tasks:
        logger.info(f"Applying task '{task}' for '{method}'")
        if task == "shapes":
            sdata_main = build_shapes(sdata, sdata_main, method, write_to_disk=False)
        if task == "cell_types":
            sdata_main = add_cell_type_annotation(
                sdata_main, sdata_path, method, write_to_disk=False
            )
        elif task == "ficture":
            sdata_main = add_ficture(
                sdata,
                sdata_main,
                method,
                ficture_arguments,
                n_ficture,
                var,
                write_to_disk=False,
            )
        elif task == "volume":
            sdata_main = calculate_volume(
                method, sdata_main, sdata_path, write_to_disk=False
            )
        elif task == "adata":
            if len(sdata.tables) == 1:
                sdata_main[f"adata_{method}"] = sdata[
                    list(sdata.tables.keys())[0]
                ].copy()
                transform_adata(sdata_main, method, data_path=data_path)

                if "cell_type_annotation" in os.listdir(
                    join(sdata_path, "results", method)
                ):
                    sdata_main = add_cell_type_annotation(
                        sdata_main, sdata_path, method, write_to_disk=False
                    )
                else:
                    logger.warning(
                        f"No cell type annotation found for '{method}'. Skipping annotation."
                    )
                if "volume" not in sdata_main[f"adata_{method}"].obs.columns:
                    sdata_main = calculate_volume(
                        method, sdata_main, sdata_path, write_to_disk=False
                    )

                if len(ficture_arguments) > 0:
                    sdata_main = add_ficture(
                        sdata,
                        sdata_main,
                        method,
                        ficture_arguments,
                        n_ficture,
                        var,
                        write_to_disk=False,
                    )
            elif len(sdata.tables) > 1:
                if logger:
                    logger.warning(
                        "Multiple adata files found for {}. Skipping adata import".format(
                            method
                        )
                    )
                else:
                    print(
                        f"Multiple table files found for {method}. Skipping adata import."
                    )
            else:
                if logger:
                    logger.warning(
                        "No adata files found for {}. Skipping adata import".format(
                            method
                        )
                    )
                else:
                    print(f"Table file missing for {method}. Skipping adata import.")
    if "shapes" in tasks:
        update_element(sdata_main, f"boundaries_{method}")
        logger.info(f"Completed shape update for '{method}'")
    elif len(set(["cell_types", "ficture", "adata", "volume"]) & tasks) > 0:
        update_element(sdata_main, f"adata_{method}")
        logger.info(f"Completed adata update for '{method}'")
