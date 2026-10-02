#!/usr/bin/env python
import argparse
import pathlib
import yaml

parser = argparse.ArgumentParser(description="Generate ComSeg sbatch scripts.")
parser.add_argument("staining", help="Staining of prior cellpose segmentation.")
parser.add_argument("CP_version", help="Cellpose version to use.")
args = parser.parse_args()
BASE_PATH = (pathlib.Path(__file__).parents[2] / "data").resolve()

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

SBATCH_DIR = BASE_PATH / f"misc/sbatches/sbatch_ComSeg_CP{args.CP_version}_{args.staining}"
SBATCH_DIR.mkdir(parents=False, exist_ok=True)

for key, value in data.items():
    cp_tag = "1" if args.staining == "nuclei" else args.CP_version
    base_segmentation = (
        f"Cellpose_1_{args.staining}_model"
        if args.staining == "nuclei"
        else f"Cellpose_{args.CP_version}_DAPI_{args.staining}"
    )

    with open(SBATCH_DIR / f"{key}.sbatch", "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 2-00:00:00
#SBATCH --mem=240G
#SBATCH --cpus-per-task=30
#SBATCH -J ComSeg_{key}_CP{cp_tag}_{args.staining}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/%x.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/%x.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new2.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
CP_VERSION="{cp_tag}"
STAINING="{args.staining}"
INPUT_PATH="{value["path"]}"
RESULT_DIR="{BASE_PATH}/samples/{key}/results/ComSeg_{base_segmentation}"
CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/comseg_sopa.py \\"${{INPUT_PATH}}\\" ${{KEY}} {base_segmentation}"
start_run_log

mamba activate segmentation
pip install -q comseg==1.8.5

mkdir -p "${{RESULT_DIR}}"
python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/comseg_sopa.py \\
  "${{INPUT_PATH}}" \\
  "${{KEY}}" \\
  "{base_segmentation}"
""")
