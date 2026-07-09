#!/usr/bin/env python
import pathlib
import yaml

BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
YAML = BASE_PATH / "misc/sample_metadata.yaml"
SBATCH_DIR = f"{BASE_PATH}/misc/sbatches/sbatch_visium"
OUT_DIR = f"{BASE_PATH}/samples/{{k}}/results/Negative_Control_Visium"

with open(YAML) as f:
    data = yaml.safe_load(f)

pathlib.Path(SBATCH_DIR).mkdir(exist_ok=True)

for k, v in data.items():
    out = OUT_DIR.format(k=k)

    text = f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 8:00:00
#SBATCH --mem=128G
#SBATCH -J visium_{k}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/visium_{k}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/visium_{k}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark.sqsh"

mamba activate segmentation
mkdir -p {out}

python $HOME/gitrepos/cellseg-benchmark/scripts/segmentation/visium_segmentation.py {v["path"]} {out}
"""

    pathlib.Path(SBATCH_DIR, f"{k}.sbatch").write_text(text)
