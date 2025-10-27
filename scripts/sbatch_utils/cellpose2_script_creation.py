import argparse
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(description="scripts for Cellpose 2 segmentation.")
parser.add_argument("staining", help="Staining of prior cellpose segmentation.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Cellpose_2"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Cellpose_2/{key}_{args.staining}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH -J CP2_{key}_{args.staining}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/CP2_{key}_{args.staining}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/CP2_{key}_{args.staining}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

mamba activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Cellpose_2_DAPI_{args.staining}
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/segmentation/cellpose_2.py {value["path"]} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Cellpose_2_DAPI_{args.staining} {args.staining}
            """)
    f.close()
