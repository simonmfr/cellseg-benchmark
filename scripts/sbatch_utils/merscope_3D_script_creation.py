import argparse
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(
    description="scripts for creating sdatas out of vpt 3D pipeline."
)
parser.add_argument("staining", help="Staining of segmentation.")
parser.add_argument("--staining2", default=None, help="Output directory.")
args = parser.parse_args()

adapt = ""
if args.staining2 is not None:
    adapt = f"_{args.staining2}"

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.save_load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_3D"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_3D/{key}_DAPI_{args.staining}{adapt}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 00:05:00
#SBATCH --mem=16G
#SBATCH -J vpt_3D_{key}_DAPI_{args.staining}{adapt}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/vpt_3D_{key}_DAPI_{args.staining}{adapt}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/vpt_3D_{key}_DAPI_{args.staining}{adapt}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

mamba activate sopa
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/merscope_3D_sdata.py {value["path"]} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/vpt_3D_DAPI_{args.staining}{adapt}
            """)
    f.close()
