#!/usr/bin/env python
import argparse
import pathlib
import yaml

parser = argparse.ArgumentParser(description="scripts for Cellpose 2 segmentation.")
parser.add_argument("staining", help="Staining of prior cellpose segmentation.")
args = parser.parse_args()
BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

pathlib.Path(
    f"{BASE_PATH}/misc/sbatches/sbatch_Cellpose_2"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"{BASE_PATH}/misc/sbatches/sbatch_Cellpose_2/{key}_{args.staining}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH -J CP2_{key}_{args.staining}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/CP2_{key}_{args.staining}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/CP2_{key}_{args.staining}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
STAINING="{args.staining}"
INPUT_PATH="{value["path"]}"

RESULT_DIR="{BASE_PATH}/samples/{key}/results/Cellpose_2_DAPI_{args.staining}"

CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/cellpose_2.py \\"${{INPUT_PATH}}\\" \\"${{RESULT_DIR}}\\" ${{STAINING}}"
start_run_log


mamba activate segmentation
mkdir -p "${{RESULT_DIR}}"
python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/cellpose_2.py \\
  "${{INPUT_PATH}}" \\
  "${{RESULT_DIR}}" \\
  "${{STAINING}}"
""")
    f.close()
