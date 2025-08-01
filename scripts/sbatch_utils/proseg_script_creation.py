import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Prepare scripts for ProSeg with prior segmentation."
)
parser.add_argument("staining", help="Staining of prior segmentation.")
parser.add_argument("CP_version", choices=["1", "2"], help="Cellpose version.")
parser.add_argument("--voxel", default=1, type=int, help="intensity ratio.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Proseg_CP{args.CP_version}_{args.staining}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    if args.staining == "nuclei":
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Proseg_CP{args.CP_version}_{args.staining}/{key}_{args.voxel}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 18:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=30
#SBATCH -J Proseg_{key}_CP1_{args.staining}_vxl_{args.voxel}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/Proseg_{key}_CP1_{args.staining}_vxl_{args.voxel}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/Proseg_{key}_CP1_{args.staining}_vxl_{args.voxel}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

source ~/.bashrc
conda activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Proseg_Cellpose_1_{args.staining}_model
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/api_proseg.py {value} {key} \
Cellpose_1_{args.staining}_model --voxel-layers {args.voxel}
            """)
        f.close()
    else:
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Proseg_CP{args.CP_version}_{args.staining}/{key}_vxl_{args.voxel}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 18:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=30
#SBATCH -J Proseg_{key}_CP{args.CP_version}_{args.staining}_vxl_{args.voxel}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/Proseg_{key}_CP{args.CP_version}_{args.staining}_vxl_{args.voxel}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/Proseg_{key}_CP{args.CP_version}_{args.staining}_vxl_{args.voxel}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

mamba activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Proseg_Cellpose_{args.CP_version}_DAPI_{args.staining}
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/api_proseg.py {value} {key} \
Cellpose_{args.CP_version}_DAPI_{args.staining} --voxel-layers {args.voxel}
            """)
        f.close()
