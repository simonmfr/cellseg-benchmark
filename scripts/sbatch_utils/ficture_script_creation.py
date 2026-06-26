#!/usr/bin/env python
import argparse
import pathlib

import yaml

BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

parser = argparse.ArgumentParser(description="Generate sbatch scripts for FICTURE.")
parser.add_argument(
    "cohort", help="Cohort name (filters samples by metadata key prefix)."
)
args = parser.parse_args()

runtime = "08:00:00" if args.cohort in ("VizgenMouseBrain", "ABCAtlas") else "04:00:00"

SBATCH_DIR = pathlib.Path(f"{BASE_PATH}/misc/sbatches/sbatch_Ficture")
SBATCH_DIR.mkdir(parents=True, exist_ok=True)

with open(f"{BASE_PATH}/misc/sample_metadata.yaml") as f:
    samples = yaml.safe_load(f)

count = 0
for key, value in samples.items():
    if not key.startswith(args.cohort):
        continue
    count += 1
    (SBATCH_DIR / f"{key}.sbatch").write_text(f"""#!/bin/bash
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t {runtime}
#SBATCH -J Ficture_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/Ficture_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/Ficture_{key}.err
#SBATCH --get-user-env
source $HOME/.bashrc
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate ficture
bash /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/ficture/ficture.sh {key} {value["path"]}/detected_transcripts.csv
""")

print(f"Wrote {count} sbatch scripts to {SBATCH_DIR}")