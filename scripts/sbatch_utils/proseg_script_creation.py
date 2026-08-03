#!/usr/bin/env python
import argparse
import pathlib
import yaml

parser = argparse.ArgumentParser(
    description="Prepare scripts for ProSeg with prior segmentation."
)
parser.add_argument("staining", help="Staining of prior segmentation.")
parser.add_argument("CP_version", choices=["1", "2"], help="Cellpose version.")
parser.add_argument("--voxel", default=1, type=int, help="intensity ratio.")
args = parser.parse_args()
BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

pathlib.Path(
    f"{BASE_PATH}/misc/sbatches/sbatch_Proseg_CP{args.CP_version}_{args.staining}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    if args.staining == "nuclei":
        f = open(
            f"{BASE_PATH}/misc/sbatches/sbatch_Proseg_CP{args.CP_version}_{args.staining}/{key}_{args.voxel}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 18:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=30
#SBATCH -J Proseg_{key}_CP1_{args.staining}_vxl_{args.voxel}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/Proseg_{key}_CP1_{args.staining}_vxl_{args.voxel}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/Proseg_{key}_CP1_{args.staining}_vxl_{args.voxel}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
CP_VERSION="1"
STAINING="{args.staining}"
VOXEL="{args.voxel}"
INPUT_PATH="{value["path"]}"

RESULT_DIR="{BASE_PATH}/samples/{key}/results/Proseg_Cellpose_1_{args.staining}_model"

CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/proseg.py \\"${{INPUT_PATH}}\\" ${{KEY}} Cellpose_1_${{STAINING}}_model --voxel-layers ${{VOXEL}}"
start_run_log

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/proseg.py \\
  "${{INPUT_PATH}}" \\
  "${{KEY}}" \\
  "Cellpose_1_${{STAINING}}_model" \\
  --voxel-layers "${{VOXEL}}"
""")
        f.close()
    else:
        f = open(
            f"{BASE_PATH}/misc/sbatches/sbatch_Proseg_CP{args.CP_version}_{args.staining}/{key}_vxl_{args.voxel}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 18:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=30
#SBATCH -J Proseg_{key}_CP{args.CP_version}_{args.staining}_vxl_{args.voxel}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/Proseg_{key}_CP{args.CP_version}_{args.staining}_vxl_{args.voxel}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/Proseg_{key}_CP{args.CP_version}_{args.staining}_vxl_{args.voxel}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
CP_VERSION="{args.CP_version}"
STAINING="{args.staining}"
VOXEL="{args.voxel}"
INPUT_PATH="{value["path"]}"

RESULT_DIR="{BASE_PATH}/samples/{key}/results/Proseg_Cellpose_{args.CP_version}_DAPI_{args.staining}"

CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/proseg.py \\"${{INPUT_PATH}}\\" ${{KEY}} Cellpose_${{CP_VERSION}}_DAPI_${{STAINING}} --voxel-layers ${{VOXEL}} --output-cell-polygon-layers cell-polygons.geojson.gz"
start_run_log

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/proseg.py \\
  "${{INPUT_PATH}}" \\
  "${{KEY}}" \\
  "Cellpose_${{CP_VERSION}}_DAPI_${{STAINING}}" \\
  --voxel-layers "${{VOXEL}}" \\
  --output-cell-polygon-layers cell-polygons.geojson.gz
""")
        f.close()
