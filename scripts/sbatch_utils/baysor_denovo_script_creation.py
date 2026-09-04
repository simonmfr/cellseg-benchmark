#!/usr/bin/env python
import argparse
import pathlib

import yaml

parser = argparse.ArgumentParser(
    description="scripts for Baysor segmentation without prior (native CLI)."
)
parser.add_argument("dimension", choices=["2D", "3D"], help="segmentation mode.")
parser.add_argument("--mem", default="90G", help="Memory per job. Serial caps at 100G per user.")
parser.add_argument("--cluster", default="serial", help="SLURM cluster.")
parser.add_argument("--partition", default="serial_long", help="SLURM partition.")
parser.add_argument("--time", default="4-00:00:00", help="Walltime per job.")
parser.add_argument("--cpus", default="8", help="Cores per job.")
args = parser.parse_args()
BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
REPO = pathlib.Path("/dss/dsshome1/0C/ra98gaq/git/cellseg-benchmark")
METHOD = f"Baysor_{args.dimension}_denovo"

with open(f"{BASE_PATH}/misc/sample_metadata.yaml") as f:
    data = yaml.safe_load(f)

pathlib.Path(f"{BASE_PATH}/misc/sbatches/sbatch_{METHOD}").mkdir(
    parents=False, exist_ok=True
)

for key, value in data.items():
    with open(f"{BASE_PATH}/misc/sbatches/sbatch_{METHOD}/{key}.sbatch", "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH --clusters={args.cluster}
#SBATCH --partition={args.partition}
#SBATCH -t {args.time}
#SBATCH --mem={args.mem}
#SBATCH --cpus-per-task={args.cpus}
#SBATCH -J {METHOD}_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/{METHOD}_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/{METHOD}_{key}.err
#SBATCH --get-user-env
#SBATCH --export=NONE

set -euo pipefail
source {REPO}/scripts/sbatch_utils/run_log.sh

KEY="{key}"
DIMENSION="{args.dimension}"
INPUT_PATH="{value["path"]}"
RESULT_DIR="{BASE_PATH}/samples/{key}/results/{METHOD}"
PARAMS="baysor=cpp-0.8.3,scale=5,n_clusters=10,mrf,no_prior"
CMD="python {REPO}/scripts/segmentation/baysor_denovo.py \\"${{INPUT_PATH}}\\" ${{KEY}} ${{DIMENSION}}"
start_run_log

export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK}}"
PYTHON="/dss/dsshome1/0C/ra98gaq/micromamba/envs/baysor_denovo/bin/python"

mkdir -p "${{RESULT_DIR}}"
"${{PYTHON}}" {REPO}/scripts/segmentation/baysor_denovo.py \\
  "${{INPUT_PATH}}" \\
  "${{KEY}}" \\
  "${{DIMENSION}}"
""")
