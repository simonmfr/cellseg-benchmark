import json
import sys
from pathlib import Path

staining = sys.argv[1]
CP_version = sys.argv[2]
confidence = sys.argv[3]

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_Baysor_CP{CP_version}_{staining}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    if staining == "nuclei":
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_Baysor_CP{CP_version}_{staining}/{key}_{confidence}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=25
#SBATCH -J Baysor_{key}_CP1_{staining}_{confidence}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/outputs/Baysor_{key}_CP1_{staining}_{confidence}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/errors/Baysor_{key}_CP1_{staining}_{confidence}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sopa.sqsh"

source ~/.bashrc
conda activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/samples/{key}/results/Baysor_2D_Cellpose_1_{staining}_model_{confidence}
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg_benchmark/scripts/api_baysor.py {value} Cellpose_1_{staining}_model {confidence} {key}
            """)
        f.close()
    else:
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sbatches/sbatch_Baysor_CP{CP_version}_{staining}/{key}_{confidence}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=25
#SBATCH -J Baysor_{key}_CP{CP_version}_{staining}_{confidence}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/outputs/Baysor_{key}_CP{CP_version}_{staining}_{confidence}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/logs/errors/Baysor_{key}_CP{CP_version}_{staining}_{confidence}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/misc/sopa.sqsh"

source ~/.bashrc
conda activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg_benchmark/samples/{key}/results/Baysor_2D_Cellpose_{CP_version}_DAPI_{staining}_{confidence}
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg_benchmark/scripts/api_baysor.py {value} Cellpose_{CP_version}_DAPI_{staining} {confidence} {key}
            """)
        f.close()
