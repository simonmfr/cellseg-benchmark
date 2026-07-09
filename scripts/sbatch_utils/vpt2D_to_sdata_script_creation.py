#!/usr/bin/env python
import argparse
import pathlib
import yaml

parser = argparse.ArgumentParser(
    description="Create SLURM jobs to convert vpt 2D outputs to sdatas."
)
parser.add_argument("staining", help="Name of staining (e.g. nuclei or PolyT).")
parser.add_argument(
    "--staining2",
    default=None,
    help="Optional nucleus staining for complex runs (e.g. nuclei).",
)
args = parser.parse_args()
BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
adapt = f"_{args.staining2}" if args.staining2 else ""

with open(f"{BASE_PATH}/misc/sample_metadata.yaml") as f:
    data = yaml.safe_load(f)

(BASE_PATH / "misc/sbatches/sbatch_vpt2D_to_sdata").mkdir(
    parents=False, exist_ok=True
)
for key, value in data.items():
    f = open(
        f"{BASE_PATH}/misc/sbatches/sbatch_vpt2D_to_sdata/{key}_{args.staining}{adapt}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 01:00:00
#SBATCH --mem=100G
#SBATCH -J vpt2D_{key}_{args.staining}{adapt}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/vpt2D_to_sdata_{key}_{args.staining}{adapt}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/vpt2D_to_sdata_{key}_{args.staining}{adapt}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark.sqsh"

mamba activate segmentation
python $HOME/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/vpt_2D_to_sdata.py {value["path"]} \
 {BASE_PATH}/samples/{key}/results/vpt_2D_DAPI_{args.staining}{adapt}
""")
    f.close()
