#!/usr/bin/env python
import pathlib
import yaml

BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
YAML = BASE_PATH / "misc/sample_metadata.yaml"
SBATCH_DIR = f"{BASE_PATH}/misc/sbatches/sbatch_voronoi"
OUT_DIR = f"{BASE_PATH}/samples/{{k}}/results/Negative_Control_Voronoi"

with open(YAML) as f:
    data = yaml.safe_load(f)

pathlib.Path(SBATCH_DIR).mkdir(exist_ok=True)

for k, v in data.items():
    out = OUT_DIR.format(k=k)

    text = f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 10:00:00
#SBATCH --mem=128G
#SBATCH -J voronoi_{k}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/voronoi_{k}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/voronoi_{k}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -euo pipefail
source $HOME/gitrepos/cellseg-benchmark/scripts/sbatch_utils/run_log.sh

KEY="{k}"
INPUT_PATH="{v["path"]}"
RESULT_DIR="{out}"

CMD="python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/voronoi_segmentation.py \\"${{INPUT_PATH}}\\" \\"${{RESULT_DIR}}\\""
start_run_log

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/voronoi_segmentation.py \\
  "${{INPUT_PATH}}" \\
  "${{RESULT_DIR}}"
"""

    pathlib.Path(SBATCH_DIR, f"{k}.sbatch").write_text(text)
