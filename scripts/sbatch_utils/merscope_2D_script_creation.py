import argparse
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(
    description="scripts for creating sdatas out of vpt 2D pipeline."
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
    data = yaml.safe_load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_2D"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_2D/{key}_DAPI_{args.staining}{adapt}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 12:00:00
#SBATCH --mem=300G
#SBATCH -J vpt_2D_{key}_DAPI_{args.staining}{adapt}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/vpt_2D_{key}_DAPI_{args.staining}{adapt}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/vpt_2D_{key}_DAPI_{args.staining}{adapt}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/benchmark.sqsh"

mamba activate segmentation
"$CONDA_PREFIX/bin/time" -v \
  -o "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/vpt_2D_{key}_DAPI_{args.staining}{adapt}_$(date +%Y%m%d_%H%M%S).time" \
python ~/gitrepos/cellseg-benchmark/scripts/segmentation/merscope_2D_sdata.py {value["path"]} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/vpt_2D_DAPI_{args.staining}{adapt}
""")
    f.close()
