#!/usr/bin/env python
import argparse
import pathlib

import yaml

parser = argparse.ArgumentParser(description="scripts for Cellpose-SAM segmentation.")
parser.add_argument(
    "cohort", nargs="?", default="", help="Cohort prefix filter. Defaults to all."
)
parser.add_argument(
    "--staining",
    help="Staining used in addition to DAPI, e.g. PolyT. Default DAPI only.",
)
args = parser.parse_args()

BASE_PATH = (pathlib.Path(__file__).parents[2] / "data").resolve()
SBATCH_DIR = BASE_PATH / "misc/sbatches/sbatch_CellposeSAM"
SBATCH_DIR.mkdir(parents=True, exist_ok=True)

METHOD = f"CellposeSAM_DAPI_{args.staining}" if args.staining else "CellposeSAM_DAPI"
STAINING_FLAG = f" --staining {args.staining}" if args.staining else ""

with open(f"{BASE_PATH}/misc/sample_metadata.yaml") as f:
    samples = yaml.safe_load(f)

count = 0
for key, value in samples.items():
    if not key.startswith(args.cohort):
        continue
    count += 1
    result_dir = f"{BASE_PATH}/samples/{key}/results/{METHOD}"
    (SBATCH_DIR / f"{key}_{METHOD}.sbatch").write_text(f"""#!/bin/bash
#SBATCH -p lrz-hgx-h100-94x4,lrz-hgx-a100-80x4,lrz-dgx-a100-80x8
#SBATCH --gres=gpu:1
#SBATCH -t 10:00:00
#SBATCH --mem=150G
#SBATCH -J {METHOD}_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/{METHOD}_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/{METHOD}_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new2.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
METHOD="{METHOD}"
INPUT_PATH="{value["path"]}"
RESULT_DIR="{result_dir}"
CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/cellpose_sam.py \\"${{INPUT_PATH}}\\" \\"${{RESULT_DIR}}\\"{STAINING_FLAG}"
start_run_log

export CELLPOSE_LOCAL_MODELS_PATH=/opt/cellpose_models
mamba activate cellposesam
mkdir -p "${{RESULT_DIR}}"

eval "${{CMD}}"
""")

print(f"Wrote {count} sbatch scripts to {SBATCH_DIR}")
