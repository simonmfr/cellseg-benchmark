from pathlib import Path

import yaml

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_transcript_tif"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_transcript_tif/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 12:00:00
#SBATCH --mem=100G
#SBATCH -J transcript_tif_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/transcript_tif_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/transcript_tif_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/benchmark_py3_12.sqsh"

cd ~/gitrepos/spatialdata
git pull -q
cd ~/gitrepos/cellseg-benchmark
git pull -q

mamba activate sopa
python scripts/transcript_tif.py {value["path"]}
""")
    f.close()
