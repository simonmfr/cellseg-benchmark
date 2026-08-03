#!/usr/bin/env python
import argparse
import pathlib
import yaml

parser = argparse.ArgumentParser(
    description="Prepare scripts for square segmentations."
)
parser.add_argument("width", type=int, help="patch width.")
parser.add_argument("unit", choices=["pixel", "microns"], help="unit of measure.")
parser.add_argument("overlap", type=int, help="patch overlap.")
parser.add_argument(
    "-ir", "--intens_rat", default=0.1, type=float, help="intensity ratio."
)
args = parser.parse_args()
BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

pathlib.Path(
    f"{BASE_PATH}/misc/sbatches/sbatch_rastered_{args.width}{args.unit}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"{BASE_PATH}/misc/sbatches/sbatch_rastered_{args.width}{args.unit}/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 08:00:00
#SBATCH --mem=128G
#SBATCH -J rastered{args.width}_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/rastered{args.width}_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/rastered{args.width}_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
INPUT_PATH="{value["path"]}"

WIDTH="{args.width}"
OVERLAP="{args.overlap}"
UNIT="{args.unit}"
INTENS_RAT="{args.intens_rat}"

RESULT_DIR="{BASE_PATH}/samples/{key}/results/Negative_Control_Rastered_{args.width}"

CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/rastered_segmentation.py \\"${{INPUT_PATH}}\\" \\"${{RESULT_DIR}}\\" ${{WIDTH}} ${{OVERLAP}} ${{UNIT}} ${{INTENS_RAT}}"
start_run_log

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/rastered_segmentation.py \\
  "${{INPUT_PATH}}" \\
  "${{RESULT_DIR}}" \\
  "${{WIDTH}}" "${{OVERLAP}}" "${{UNIT}}" "${{INTENS_RAT}}"
""")
    f.close()
