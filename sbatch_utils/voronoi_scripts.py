import json
from pathlib import Path

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatch_voronoi").mkdir(
    parents=False, exist_ok=True
)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatch_voronoi/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
            
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -J voronoi_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/voronoi_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/voronoi_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

mamba activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Negative_Control_Voronoi
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/voronoi_segmentation.py \
 {value} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Negative_Control_Voronoi 100000
            """)
    f.close()
