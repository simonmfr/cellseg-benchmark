#!/usr/bin/env python
import pathlib
import yaml

BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

pathlib.Path(
    f"{BASE_PATH}/misc/sbatches/sbatch_nuclei_model"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"{BASE_PATH}/misc/sbatches/sbatch_nuclei_model/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH -J nuclei_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/nuclei_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/nuclei_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{key}"
INPUT_PATH="{value["path"]}"

RESULT_DIR="{BASE_PATH}/samples/{key}/results/Cellpose_1_nuclei_model"

CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/nuclei.py \\"${{INPUT_PATH}}\\" \\"${{RESULT_DIR}}\\""
start_run_log

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/nuclei.py \\
  "${{INPUT_PATH}}" \\
  "${{RESULT_DIR}}"
""")
    f.close()
