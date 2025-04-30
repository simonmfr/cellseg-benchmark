import json
import sys
from pathlib import Path

staining = sys.argv[1]
CP_version = sys.argv[2]

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_ComSeg_CP{CP_version}_{staining}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    if staining == "nuclei":
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_ComSeg_CP{CP_version}_{staining}/{key}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 2-00:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=10
#SBATCH -J ComSeg_{key}_CP1_{staining}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/outputs/ComSeg_{key}_CP1_{staining}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/errors/ComSeg_{key}_CP1_{staining}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sopa.sqsh"

source ~/.bashrc
conda activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/samples/{key}/results/Baysor_2D_Cellpose_1_{staining}_model
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg_benchmark/scripts/api_comseg.py {value} {key} \
Cellpose_1_{staining}_model
            """)
        f.close()
    else:
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_ComSeg_CP{CP_version}_{staining}/{key}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 2-00:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=10
#SBATCH -J ComSeg_{key}_CP{CP_version}_{staining}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/outputs/ComSeg_{key}_CP{CP_version}_{staining}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/errors/ComSeg_{key}_CP{CP_version}_{staining}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sopa.sqsh"

mamba activate sopa
pip install igraph
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/samples/{key}/results/ComSeg_Cellpose_{CP_version}_DAPI_{staining}
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg_benchmark/scripts/api_comseg.py {value} {key} \
 Cellpose_{CP_version}_DAPI_{staining}
            """)
        f.close()
