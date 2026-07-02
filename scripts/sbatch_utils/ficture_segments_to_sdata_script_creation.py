#!/usr/bin/env python
"""Create one SLURM job per sample to build FICTURE cell boundaries as sdatas.

Each job runs scripts/ficture/ficture_segments_to_sdata.py, which writes the
raw FICTURE segments (Ficture_segments) and the nuclei-split cells
(Ficture_segments_dapi) as SpatialData .zarr objects, and aggregates MERSCOPE
transcripts into a cell x gene table on the split cells.
"""
import argparse
import pathlib

import yaml

BASE_PATH = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
CPUS = 16  # per-factor parallelism (exposed to the script as SLURM_CPUS_PER_TASK)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--res", type=float, default=1.5, help="grid size in um")
args = parser.parse_args()

with open(f"{BASE_PATH}/misc/sample_metadata.yaml") as f:
    data = yaml.safe_load(f)

out_dir = pathlib.Path(f"{BASE_PATH}/misc/sbatches/sbatch_ficture_segments_to_sdata")
out_dir.mkdir(parents=False, exist_ok=True)

for key, value in data.items():
    with open(out_dir / f"{key}.sbatch", "w") as f:
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task={CPUS}
#SBATCH --mem=128G
#SBATCH -J ficture_segments_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/ficture_segments_to_sdata_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/ficture_segments_to_sdata_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark.sqsh"

mamba activate segmentation
python ~/gitrepos/cellseg-benchmark/scripts/ficture/ficture_segments_to_sdata.py \\
 {key} --res {args.res} --data-path {value["path"]}
""")
