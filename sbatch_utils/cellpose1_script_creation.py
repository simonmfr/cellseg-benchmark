import json
from pathlib import Path
import sys

staining = sys.argv[1]

with open("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json") as f:
    data = json.load(f)

Path(f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatch_Cellpose_1").mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatch_Cellpose_1/{key}_{staining}.sbatch", "w")
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=30
#SBATCH -J CP1_{key}_{staining}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/CP1_{key}_{staining}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/CP1_{key}_{staining}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

source ~/.bashrc
conda activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Cellpose_1_DAPI_{staining}
python /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa_api/api_cellpose_1.py {value} /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Cellpose_1_DAPI_{staining} {staining}
            """)
    f.close()
