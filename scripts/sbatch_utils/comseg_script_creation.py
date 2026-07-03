#!/usr/bin/env python
import argparse
import pathlib
import yaml

parser = argparse.ArgumentParser(description="Generate ComSeg sbatch scripts.")
parser.add_argument("staining", help="Staining of prior cellpose segmentation.")
parser.add_argument("CP_version", help="Cellpose version to use.")
args = parser.parse_args()
BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

outdir = pathlib.Path(
    f"{BASE_PATH}/misc/sbatches/"
    f"sbatch_ComSeg_CP{args.CP_version}_{args.staining}"
)
outdir.mkdir(parents=False, exist_ok=True)

for key, value in data.items():
    cp_tag = "CP1" if args.staining == "nuclei" else f"CP{args.CP_version}"
    job_name = f"ComSeg_{key}_{cp_tag}_{args.staining}"
    result_dir = (
        f"Cellpose_{args.CP_version}_DAPI_{args.staining}"
        if args.staining != "nuclei"
        else f"Cellpose_1_{args.staining}_model"
    )

    sbatch_path = outdir / f"{key}.sbatch"
    with open(sbatch_path, "w") as f:
        f.write(
            f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 12:00:00
#SBATCH --mem=160G
#SBATCH --cpus-per-task=20
#SBATCH -J {job_name}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/%x.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/%x.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_py3_12.sqsh"

set -eu

cd $HOME/gitrepos/spatialdata
git pull -q
cd $HOME/gitrepos/cellseg-benchmark
git pull -q

set -eu

mamba activate sopa

mkdir -p {BASE_PATH}/samples/{key}/results/ComSeg_{result_dir}

python scripts/segmentation/comseg.py {value["path"]} {key} {result_dir}
"""
        )
