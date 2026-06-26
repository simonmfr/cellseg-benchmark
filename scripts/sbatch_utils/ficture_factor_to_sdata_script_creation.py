#!/usr/bin/env python
import argparse
import pathlib
import yaml

BASE_PATH = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"

parser = argparse.ArgumentParser(
    description="Create SLURM jobs to convert Ficture factor x gene counts to sdatas."
)
args = parser.parse_args()

with open(f"{BASE_PATH}/misc/sample_metadata.yaml") as f:
    data = yaml.safe_load(f)

pathlib.Path(f"{BASE_PATH}/misc/sbatches/sbatch_ficture_factor_to_sdata").mkdir(
    parents=False, exist_ok=True
)
for key, value in data.items():
    out = f"{BASE_PATH}/samples/{key}/results/Ficture"
    f = open(
        f"{BASE_PATH}/misc/sbatches/sbatch_ficture_factor_to_sdata/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 00:20:00
#SBATCH --mem=16G
#SBATCH -J ficture_factor_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/ficture_factor_to_sdata_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/ficture_factor_to_sdata_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark.sqsh"

mamba activate segmentation
python ~/gitrepos/cellseg-benchmark/scripts/ficture/ficture_factor_to_h5ad.py \
 {out}/output/model.posterior.count.tsv.gz \
 --zarr {out}/../Ficture_factors/sdata.zarr
""")
    f.close()