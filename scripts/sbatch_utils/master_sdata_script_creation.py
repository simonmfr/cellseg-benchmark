import json
from pathlib import Path

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_master_sdata"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_master_sdata/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 01:00:00
#SBATCH --mem=150G
#SBATCH -J master_sdata_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/master_sdata_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/master_sdata_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/cellseg_benchmark.sqsh"


cd ~/gitrepos/cellseg_benchmark
git pull
mamba activate cellseg_benchmark
python scripts/master_sdata.py {key} {value} z3 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark 
            """)
    f.close()
