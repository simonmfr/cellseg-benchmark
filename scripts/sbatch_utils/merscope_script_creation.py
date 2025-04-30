import json
from pathlib import Path

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_merscope"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_merscope/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 12:00:00
#SBATCH --mem=300G
#SBATCH -J merscope_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/outputs/merscope_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/errors/merscope_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sopa.sqsh"

mamba activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/samples/{key}/results/Cellpose_1_Merlin
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg_benchmark/scripts/merscope_sdata.py {value} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/samples/{key}/results/Cellpose_1_Merlin
            """)
    f.close()
