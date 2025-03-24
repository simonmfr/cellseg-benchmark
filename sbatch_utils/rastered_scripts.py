import json
import sys
from pathlib import Path

size = sys.argv[1]
unit = sys.argv[2]
overlap = sys.argv[3]
if 1 >= int(sys.argv[4]) >= 0:
    intens_rat = int(sys.argv[4])
else:
    intens_rat = 0.1

assert unit=="microns" or unit=="pixel", "unit must be either microns or pixel."

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatch_rastered_{size}{unit}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatch_rastered_{size}{unit}/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
            
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128G
#SBATCH -J rastered{size}_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/rastered{size}_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/rastered{size}_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

mamba activate sopa
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/rastered_segmentation.py \
 {value} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Negative_Control_Rastered_{size} \
 {size} {overlap} {unit} {intens_rat}
            """)
    f.close()
