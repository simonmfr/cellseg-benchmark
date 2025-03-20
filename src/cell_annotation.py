import os
import sys
import warnings
from subprocess import run

import scanpy as sc
from sopa.io.explorer import write_cell_categories

import cell_annotation_utils

warnings.filterwarnings("ignore")

sample_name = sys.argv[1]
method_name = sys.argv[2]
data_dir = sys.argv[3]
mad_factor = int(sys.argv[4]) if int(sys.argv[4]) > 0 else 3

# Setup and unpacking
if "SLURM_CPUS_PER_TASK" in os.environ:
    sc.settings.n_jobs = int(os.environ["SLURM_CPUS_PER_TASK"])
    print(sc.settings.n_jobs)

path = os.path.join(data_dir, "samples", sample_name, "results", method_name)
allen_mmc_dir = os.path.join(path, "cell_type_annotation")
mmc_out = [
    f
    for f in os.listdir(allen_mmc_dir)
    if "_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping" in f
]
if len(mmc_out) > 1:
    AssertionError("More than one hierarchical Mapping file found")
json_name = os.path.splitext(mmc_out[0])[0]

for file in mmc_out:
    if file.endswith(".zip"):
        file_path = os.path.join(allen_mmc_dir, file)
        unzip_path = os.path.join(allen_mmc_dir, file.replace(".zip", ""))
        run(["unzip", file_path, "-d", unzip_path])
        run(["rm", file_path])

adata = cell_annotation_utils.celltype_mapping(
    path, json_name, allen_mmc_dir, data_dir, mad_factor
)

write_cell_categories(os.path.join(path, "sdata.explorer"), adata)
