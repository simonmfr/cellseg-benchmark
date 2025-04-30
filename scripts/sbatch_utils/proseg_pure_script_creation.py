import json
import sys
from pathlib import Path

voxel = sys.argv[1]

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_Proseg_pure"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_Proseg_pure/{key}_vxl_{voxel}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 18:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=30
#SBATCH -J Proseg_pure_{key}_vxl_{voxel}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/outputs/Proseg_pure_{key}_vxl_{voxel}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/errors/Proseg_pure_{key}_vxl_{voxel}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sopa.sqsh"

mamba activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/samples/{key}/results/Proseg_pure
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg_benchmark/scripts/proseg_pure.py {value} {key} --voxel-layers {voxel}
            """)
    f.close()
