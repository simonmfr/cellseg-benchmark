import json
import sys
from pathlib import Path

staining = sys.argv[1]
CP_version = sys.argv[2]
voxel = sys.argv[3]

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_Proseg_CP{CP_version}_{staining}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    if staining == "nuclei":
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_Proseg_CP{CP_version}_{staining}/{key}_{voxel}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 18:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=30
#SBATCH -J Proseg_{key}_CP1_{staining}_vxl_{voxel}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/outputs/Proseg_{key}_CP1_{staining}_vxl_{voxel}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/errors/Proseg_{key}_CP1_{staining}_vxl_{voxel}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sopa.sqsh"

source ~/.bashrc
conda activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/samples/{key}/results/Proseg_Cellpose_1_{staining}_model
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg_benchmark/scripts/api_proseg.py {value} {key} \
Cellpose_1_{staining}_model  --voxel-layers {voxel}
            """)
        f.close()
    else:
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_Proseg_CP{CP_version}_{staining}/{key}_vxl_{voxel}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 18:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=30
#SBATCH -J Proseg_{key}_CP{CP_version}_{staining}_vxl_{voxel}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/outputs/Proseg_{key}_CP{CP_version}_{staining}_vxl_{voxel}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/errors/Proseg_{key}_CP{CP_version}_{staining}_vxl_{voxel}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sopa.sqsh"

mamba activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/samples/{key}/results/Proseg_Cellpose_{CP_version}_DAPI_{staining}
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg_benchmark/scripts/api_proseg.py {value} {key} \
Cellpose_{CP_version}_DAPI_{staining}  --voxel-layers {voxel}
            """)
        f.close()
