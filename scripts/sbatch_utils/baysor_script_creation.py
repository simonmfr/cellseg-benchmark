import argparse
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(description="scripts for Cellpose 2 segmentation.")
parser.add_argument("staining", help="Staining of prior cellpose segmentation.")
parser.add_argument("CP_version", help="Cellpose version of prior  segmentation.")
parser.add_argument("confidence", help="Confidence of prior cellpose segmentation.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Baysor_CP{args.CP_version}_{args.staining}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    if args.staining == "nuclei":
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Baysor_CP{args.CP_version}_{args.staining}/{key}_{args.confidence}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=25
#SBATCH -J Baysor_{key}_CP1_{args.staining}_{args.confidence}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/Baysor_{key}_CP1_{args.staining}_{args.confidence}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/Baysor_{key}_CP1_{args.staining}_{args.confidence}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

source ~/.bashrc
conda activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Baysor_2D_Cellpose_1_{args.staining}_model_{args.confidence}
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/segmentation/baysor.py {value["path"]} Cellpose_1_{args.staining}_model {args.confidence} {key}
""")
        f.close()
    else:
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Baysor_CP{args.CP_version}_{args.staining}/{key}_{args.confidence}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=25
#SBATCH -J Baysor_{key}_CP{args.CP_version}_{args.staining}_{args.confidence}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/Baysor_{key}_CP{args.CP_version}_{args.staining}_{args.confidence}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/Baysor_{key}_CP{args.CP_version}_{args.staining}_{args.confidence}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

source ~/.bashrc
conda activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Baysor_2D_Cellpose_{args.CP_version}_DAPI_{args.staining}_{args.confidence}
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/segmentation/baysor.py {value["path"]} Cellpose_{args.CP_version}_DAPI_{args.staining} {args.confidence} {key}
""")
        f.close()
