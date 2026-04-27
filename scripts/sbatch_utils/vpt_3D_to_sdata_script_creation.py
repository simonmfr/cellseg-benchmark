#!/usr/bin/env python

import argparse
import pathlib
import yaml
import cellseg_benchmark as cb

parser = argparse.ArgumentParser(
    description="Create SLURM jobs to convert vpt 3D outputs to sdatas."
)
parser.add_argument("staining", help="Staining of segmentation.")
parser.add_argument("--staining2", default=None)
args = parser.parse_args()

adapt = f"_{args.staining2}" if args.staining2 else ""

with open(f"{cb.BASE_PATH}/misc/sample_metadata.yaml") as f:
    data = yaml.safe_load(f)

pathlib.Path(f"{cb.BASE_PATH}/misc/sbatches/sbatch_vpt_3D_to_sdata").mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"{cb.BASE_PATH}/misc/sbatches/sbatch_vpt_3D_to_sdata/{key}_DAPI_{args.staining}{adapt}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 00:05:00
#SBATCH --mem=75G
#SBATCH -J vpt_3D_{key}_DAPI_{args.staining}{adapt}
#SBATCH -o {cb.BASE_PATH}/misc/logs/outputs/vpt_3D_to_sdata_{key}_DAPI_{args.staining}{adapt}.out
#SBATCH -e {cb.BASE_PATH}/misc/logs/errors/vpt_3D_to_sdata_{key}_DAPI_{args.staining}{adapt}.err
#SBATCH --container-image="{cb.BASE_PATH}/misc/enroot_images/benchmark.sqsh"

mamba activate segmentation
python ~/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/vpt_3D_to_sdata.py {value["path"]} \
 {cb.BASE_PATH}/samples/{key}/results/vpt_3D_DAPI_{args.staining}{adapt}
""")
    f.close()