#!/usr/bin/env python
import pathlib
import yaml

BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

pathlib.Path(
    f"{BASE_PATH}/misc/sbatches/sbatch_transcript_tif"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"{BASE_PATH}/misc/sbatches/sbatch_transcript_tif/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 12:00:00
#SBATCH --mem=100G
#SBATCH -J transcript_tif_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/transcript_tif_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/transcript_tif_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -euo pipefail

INPUT_PATH="{value["path"]}"

mamba activate segmentation

python $HOME/gitrepos/cellseg-benchmark/scripts/transcript_tif.py \\
  "${{INPUT_PATH}}"
""")
    f.close()
