#!/usr/bin/env python
import argparse
import pathlib

import yaml

parser = argparse.ArgumentParser(
    description="Create SLURM jobs to convert Baysor de novo outputs to sdatas."
)
parser.add_argument("dimension", choices=["2D", "3D"], help="segmentation mode.")
args = parser.parse_args()
BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
# DSS checkout, shared by both clusters. Switch to $HOME/gitrepos once merged.
REPO = pathlib.Path("/dss/dsshome1/0C/ra98gaq/git/cellseg-benchmark")
METHOD = f"Baysor_{args.dimension}_denovo"

with open(f"{BASE_PATH}/misc/sample_metadata.yaml") as f:
    data = yaml.safe_load(f)

(BASE_PATH / f"misc/sbatches/sbatch_{METHOD}_to_sdata").mkdir(
    parents=False, exist_ok=True
)
for key, value in data.items():
    f = open(
        f"{BASE_PATH}/misc/sbatches/sbatch_{METHOD}_to_sdata/{key}.sbatch", "w"
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 01:00:00
#SBATCH --mem=150G
#SBATCH -J {METHOD}_to_sdata_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/{METHOD}_to_sdata_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/{METHOD}_to_sdata_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

mamba activate segmentation
python {REPO}/scripts/seg_postprocessing/baysor_denovo_to_sdata.py {value["path"]} \
 {BASE_PATH}/samples/{key}/results/{METHOD}
""")
    f.close()
