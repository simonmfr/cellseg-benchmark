import json
from pathlib import Path
from sys import argv

staining = argv[1]
adapt = ""
if len(argv) >= 3:
    staining2 = argv[2]
    adapt = f"_{staining2}"

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_3D"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_3D/{key}_DAPI_{staining}{adapt}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 12:00:00
#SBATCH --mem=300G
#SBATCH -J vtp_3D_{key}_DAPI_{staining}{adapt}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/vpt_3D_{key}_DAPI_{staining}{adapt}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/vpt_3D_{key}_DAPI_{staining}{adapt}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

mamba activate sopa
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/merscope_3D_sdata.py {value} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/vpt_3D_DAPI_{staining}{adapt}
            """)
    f.close()
