#!/usr/bin/env python
import argparse
import pathlib

import yaml

parser = argparse.ArgumentParser(
    description="scripts for ProSeg 3D seeded by Baysor 3D (3DMERFISH thick sections)."
)
parser.add_argument("--mem", default="300G", help="Memory per job.")
parser.add_argument("--time", default="48:00:00", help="Walltime per job.")
parser.add_argument("--cpus", default="32", help="Cores per job.")
parser.add_argument(
    "--repo",
    default="$HOME/gitrepos/cellseg-benchmark",
    help="Repo with proseg_3D_baysor.py (default: container main).",
)
args = parser.parse_args()

BASE_PATH = (pathlib.Path(__file__).parents[2] / "data").resolve()
METHOD = "Proseg_3D_Baysor_3D_denovo"

THICKNESS = {
    "3DMERFISH_cortex": 100,
    "3DMERFISH_hypothalamus": 200,
}  # µm, one layer per 1 µm plane

with open(BASE_PATH / "misc/sample_metadata.yaml") as f:
    data = yaml.safe_load(f)

(BASE_PATH / f"misc/sbatches/sbatch_{METHOD}").mkdir(parents=False, exist_ok=True)

for key, value in data.items():
    if key not in THICKNESS:
        continue
    sigma_z = 1 / THICKNESS[key]  # 1 µm, proseg z is rescaled to [0, 1]
    with open(BASE_PATH / f"misc/sbatches/sbatch_{METHOD}/{key}.sbatch", "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t {args.time}
#SBATCH --mem={args.mem}
#SBATCH --cpus-per-task={args.cpus}
#SBATCH -J {METHOD}_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/{METHOD}_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/{METHOD}_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -euo pipefail
source "$HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh"

KEY="{key}"
INPUT_PATH="{BASE_PATH}/samples/{key}/results/Baysor_3D_denovo"
RESULT_DIR="{BASE_PATH}/samples/{key}/results/{METHOD}"
FLAGS="--voxel-layers {THICKNESS[key]} --voxel-size 2 --burnin-voxel-size 4 --diffusion-sigma-z {sigma_z:g} --diffusion-proposal-sigma-z {2 * sigma_z:g} --output-cell-polygon-layers cell-polygons-layers.geojson.gz"
PARAMS="proseg=3.0.7,prior=Baysor_3D_denovo,${{FLAGS// /,}}"
CMD="python {args.repo}/scripts/segmentation/proseg_3D_baysor.py \\"${{INPUT_PATH}}\\" \\"${{RESULT_DIR}}\\" ${{FLAGS}}"
start_run_log

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python "{args.repo}/scripts/segmentation/proseg_3D_baysor.py" \\
  "${{INPUT_PATH}}" \\
  "${{RESULT_DIR}}" \\
  ${{FLAGS}}
""")
