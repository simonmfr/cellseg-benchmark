#!/usr/bin/env python
import argparse
import pathlib

import yaml

parser = argparse.ArgumentParser(description="scripts for CellSAM segmentation.")
parser.add_argument(
    "cohort", nargs="?", default="", help="Cohort prefix filter. Defaults to all."
)
parser.add_argument(
    "--use_polyt", action="store_true", help="Use DAPI + PolyT instead of DAPI only."
)
parser.add_argument("--block_size", type=int, default=512, help="CellSAM block size.")
parser.add_argument("--patch_width", type=int, default=8192, help="Sopa patch width.")
args = parser.parse_args()

BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
SBATCH_DIR = BASE_PATH / "misc/sbatches/sbatch_CellSAM"
SBATCH_DIR.mkdir(parents=True, exist_ok=True)

METHOD = "CellSAM_DAPI_PolyT" if args.use_polyt else "CellSAM_DAPI"
POLYT_FLAG = "--use_polyt " if args.use_polyt else ""

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
#SBATCH -t 1-12:00:00
#SBATCH --mem=300G
#SBATCH -J CellSAM_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/CellSAM_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/CellSAM_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
METHOD="{METHOD}"
INPUT_PATH="{value["path"]}"
RESULT_DIR="{result_dir}"
CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/cellsam.py \\"${{INPUT_PATH}}\\" \\"${{RESULT_DIR}}\\" {POLYT_FLAG}--block_size {args.block_size} --patch_width {args.patch_width}"
start_run_log

mamba activate cellsam
mkdir -p "${{RESULT_DIR}}"

eval "${{CMD}}"
""")

print(f"Wrote {count} sbatch scripts to {SBATCH_DIR}")
