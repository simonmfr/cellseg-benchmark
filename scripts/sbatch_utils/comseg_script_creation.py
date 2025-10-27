import argparse
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(description="scripts for ComSeg segmentation.")
parser.add_argument("staining", help="Staining of prior cellpose segmentation.")
parser.add_argument("CP_version", help="Cellpose version to use.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_ComSeg_CP{args.CP_version}_{args.staining}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    if args.staining == "nuclei":
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_ComSeg_CP{args.CP_version}_{args.staining}/{key}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 2-00:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=10
#SBATCH -J ComSeg_{key}_CP1_{args.staining}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/ComSeg_{key}_CP1_{args.staining}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/ComSeg_{key}_CP1_{args.staining}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

source ~/.bashrc
conda activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Baysor_2D_Cellpose_1_{args.staining}_model
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/segmentation/comseg.py {value["path"]} {key} \
Cellpose_1_{args.staining}_model
""")
        f.close()
    else:
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_ComSeg_CP{args.CP_version}_{args.staining}/{key}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 2-00:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=10
#SBATCH -J ComSeg_{key}_CP{args.CP_version}_{args.staining}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/ComSeg_{key}_CP{args.CP_version}_{args.staining}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/ComSeg_{key}_CP{args.CP_version}_{args.staining}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

mamba activate sopa
pip install igraph
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/ComSeg_Cellpose_{args.CP_version}_DAPI_{args.staining}
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/segmentation/comseg.py {value["path"]} {key} \
 Cellpose_{args.CP_version}_DAPI_{args.staining}
""")
        f.close()
