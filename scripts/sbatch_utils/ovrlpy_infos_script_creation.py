#!/usr/bin/env python
import argparse
import pathlib
import yaml

parser = argparse.ArgumentParser(description="scripts for computing ovrlpy statistics.")
parser.add_argument(
    "--recompute", action="store_true", help="Consider genotype differentiation"
)
args = parser.parse_args()
BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

pathlib.Path(
    f"{BASE_PATH}/misc/sbatches/sbatch_ovrlpy_stats"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"{BASE_PATH}/misc/sbatches/sbatch_ovrlpy_stats/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH -q cpu
#SBATCH -t 01:30:00
#SBATCH --mem=50G
#SBATCH -J ovrlpy_stats_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/ovrlpy_stats_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/ovrlpy_stats_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark.sqsh"

mamba activate seg_postprocessing
python $HOME/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/ovrlpy_infos.py {key} {value["path"]} {"--recompute" if args.recompute else ""}
""")
    f.close()
