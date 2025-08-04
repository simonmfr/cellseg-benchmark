import json
from pathlib import Path

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Ficture"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Ficture/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=112
#SBATCH -t 02:00:00
#SBATCH -J Ficture_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/Ficture_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/Ficture_{key}.err
#SBATCH --get-user-env

source $HOME/.bashrc
conda activate ficture
bash /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/ficture.sh {key} {value}/detected_transcripts.csv
            """)
    f.close()
