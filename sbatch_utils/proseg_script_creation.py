import json
import sys
from pathlib import Path

staining = sys.argv[1]
CP_version = sys.argv[2]

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatch_Baysor_CP{CP_version}_{staining}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatch_Proseg_CP{CP_version}_{staining}/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-12:00:00
#SBATCH --mem=300G
#SBATCH -J Proseg_{key}_CP{CP_version}_{staining}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/Proseg_{key}_CP{CP_version}_{staining}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/Proseg_{key}_CP{CP_version}_{staining}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

source ~/.bashrc
conda activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Proseg_CP{CP_version}_DAPI_{staining}
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/api_proseg.py {value} {key} Cellpose_{CP_version}_DAPI_{staining}
            """)
    f.close()
