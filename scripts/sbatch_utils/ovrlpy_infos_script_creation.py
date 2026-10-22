#!/usr/bin/env python
import argparse
import pathlib
import yaml

parser = argparse.ArgumentParser(description="scripts for computing ovrlpy statistics.")
parser.add_argument("--linux_cluster", action="store_true", help="Adapt sbatch call for linux cluster.")
parser.add_argument(
    "--recompute", action="store_true", help="Consider genotype differentiation"
)
args = parser.parse_args()
BASE_PATH = (pathlib.Path(__file__).parents[2] / "data").resolve()

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

pathlib.Path(
    f"{BASE_PATH}/misc/sbatches/sbatch_ovrlpy_stats"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    if args.linux_cluster:
        f = open(
            f"{BASE_PATH}/misc/sbatches/sbatch_ovrlpy_stats_linux/{key}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 01:30:00
#SBATCH --mem=50G
#SBATCH -J ovrlpy_stats_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/ovrlpy_stats_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/ovrlpy_stats_{key}.err
#SBATCH --get-user-env

source $HOME/.bashrc
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate cellseg_benchmark

python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/seg_postprocessing/ovrlpy_infos.py {key} {value["path"]} {"--recompute" if args.recompute else ""}
""")
        f.close()
    else:
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
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

mamba activate seg_postprocessing
python $HOME/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/ovrlpy_infos.py {key} {value["path"]} {"--recompute" if args.recompute else ""}
""")
        f.close()
