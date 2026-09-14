#!/usr/bin/env python
import argparse
import pathlib
import yaml

parser = argparse.ArgumentParser(description="scripts for Baysor segmentation.")
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
    f"{BASE_PATH}/misc/sbatches/sbatch_Baysor_CP{args.CP_version}_{args.staining}"
).mkdir(parents=False, exist_ok=True)

for key, value in data.items():
    if args.staining == "nuclei":
        cp_tag = "1"
        result_dir_name = f"Baysor_2D_Cellpose_1_{args.staining}_model_{args.confidence}"
        model_token_log = "Cellpose_1_${STAINING}_model"
        model_token_cmd = '"Cellpose_1_${STAINING}"_model'
    else:
        cp_tag = args.CP_version
        result_dir_name = f"Baysor_2D_Cellpose_{args.CP_version}_DAPI_{args.staining}_{args.confidence}"
        model_token_log = "Cellpose_${CP_VERSION}_DAPI_${STAINING}"
        model_token_cmd = '"Cellpose_${CP_VERSION}_DAPI_${STAINING}"'

    with open(
        f"{BASE_PATH}/misc/sbatches/sbatch_Baysor_CP{args.CP_version}_{args.staining}/{key}_{args.confidence}.sbatch",
        "w",
    ) as f:
        f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=75G
#SBATCH -J Baysor_{key}_CP{cp_tag}_{args.staining}_{args.confidence}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/Baysor_{key}_CP{cp_tag}_{args.staining}_{args.confidence}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/Baysor_{key}_CP{cp_tag}_{args.staining}_{args.confidence}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
CP_VERSION="{args.CP_version}"
STAINING="{args.staining}"
CONFIDENCE="{args.confidence}"
INPUT_PATH="{value["path"]}"
RESULT_DIR="{BASE_PATH}/samples/{key}/results/{result_dir_name}"
CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/baysor_sopa.py \\"${{INPUT_PATH}}\\" {model_token_log} ${{CONFIDENCE}} ${{KEY}}"
start_run_log

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/baysor_sopa.py \
  "${{INPUT_PATH}}" \
  {model_token_cmd} \
  "${{CONFIDENCE}}" \
  "${{KEY}}"
""")
