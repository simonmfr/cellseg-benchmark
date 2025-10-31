import argparse
from pathlib import Path
import yaml

parser = argparse.ArgumentParser(description="Generate ComSeg sbatch scripts.")
parser.add_argument("staining", help="Staining of prior cellpose segmentation.")
parser.add_argument("CP_version", help="Cellpose version to use.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

outdir = Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/"
    f"sbatch_ComSeg_CP{args.CP_version}_{args.staining}"
)
outdir.mkdir(parents=False, exist_ok=True)

for key, value in data.items():
    cp_tag = "CP1" if args.staining == "nuclei" else f"CP{args.CP_version}"
    job_name = f"ComSeg_{key}_{cp_tag}_{args.staining}"
    result_dir = f"Cellpose_{args.CP_version}_DAPI_{args.staining}" if args.staining != "nuclei" else f"Cellpose_1_{args.staining}_model"

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
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/%x.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/%x.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

set -eu

source ~/.bashrc
conda activate sopa

mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/ComSeg_{result_dir}

python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/segmentation/comseg.py \\
    {value["path"]} \\
    {key} \\
    {result_dir}
"""
        )
