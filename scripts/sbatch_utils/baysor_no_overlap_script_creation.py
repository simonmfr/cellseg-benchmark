#!/usr/bin/env python
import argparse
import pathlib
import yaml

parser = argparse.ArgumentParser(
    description="scripts for Baysor segmentation without overlaps."
)
parser.add_argument("staining", help="Staining of prior cellpose segmentation.")
parser.add_argument("CP_version", help="Cellpose version of prior  segmentation.")
parser.add_argument("confidence", help="Confidence of prior cellpose segmentation.")
args = parser.parse_args()
BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

pathlib.Path(
    f"{BASE_PATH}/misc/sbatches/sbatch_Baysor_no_overlap_CP{args.CP_version}_{args.staining}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    if args.staining == "nuclei":
        f = open(
            f"{BASE_PATH}/misc/sbatches/sbatch_Baysor_no_overlap_CP{args.CP_version}_{args.staining}/{key}_{args.confidence}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH -J Baysor_no_overlap_{key}_CP1_{args.staining}_{args.confidence}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/Baysor_no_overlap_{key}_CP1_{args.staining}_{args.confidence}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/Baysor_no_overlap_{key}_CP1_{args.staining}_{args.confidence}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
CP_VERSION="{args.CP_version}"
STAINING="{args.staining}"
CONFIDENCE="{args.confidence}"
INPUT_PATH="{value["path"]}"

ZARR_DIR="{BASE_PATH}/samples/{key}/results/Baysor_2D_Cellpose_1_{args.staining}_model_{args.confidence}"
RESULT_DIR="{BASE_PATH}/samples/{key}/results/Baysor_2D_no_overlap_Cellpose_1_{args.staining}_model_{args.confidence}"

CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/baysor_no_overlap.py \\"${{INPUT_PATH}}\\" ${{ZARR_DIR}} ${{RESULT_DIR}}"
start_run_log

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/baysor_no_overlap.py \
  "${{INPUT_PATH}}" \
  "${{ZARR_DIR}}" \
  "${{RESULT_DIR}}"
""")
        f.close()
    else:
        f = open(
            f"{BASE_PATH}/misc/sbatches/sbatch_Baysor_no_overlap_CP{args.CP_version}_{args.staining}/{key}_{args.confidence}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=25
#SBATCH -J Baysor_no_overlap_{key}_CP{args.CP_version}_{args.staining}_{args.confidence}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/Baysor_no_overlap_{key}_CP{args.CP_version}_{args.staining}_{args.confidence}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/Baysor_no_overlap_{key}_CP{args.CP_version}_{args.staining}_{args.confidence}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
CP_VERSION="{args.CP_version}"
STAINING="{args.staining}"
CONFIDENCE="{args.confidence}"
INPUT_PATH="{value["path"]}"

ZARR_DIR="{BASE_PATH}/samples/{key}/results/Baysor_2D_Cellpose_1_{args.staining}_model_{args.confidence}"
RESULT_DIR="{BASE_PATH}/samples/{key}/results/Baysor_2D_no_overlap_Cellpose_{args.CP_version}_DAPI_{args.staining}_{args.confidence}"

CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/baysor_no_overlap.py \\"${{INPUT_PATH}}\\" ${{ZARR_DIR}} ${{RESULT_DIR}}"
start_run_log

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/baysor_no_overlap.py "${{INPUT_PATH}}" "${{ZARR_DIR}}" "${{RESULT_DIR}}"
""")
        f.close()
