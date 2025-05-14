import json
from pathlib import Path

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_ficture_stats"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_ficture_stats/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 03:00:00
#SBATCH --mem=200G
#SBATCH -J ficture_stats_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/ficture_stats_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/ficture_stats_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/cellseg_benchmark.sqsh"

cd ~/gitrepos/cellseg-benchmark
git pull
checkout Ficture
mamba activate cellseg_benchmark
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/ficture_infos.py {key} {value}
""")
    f.close()
