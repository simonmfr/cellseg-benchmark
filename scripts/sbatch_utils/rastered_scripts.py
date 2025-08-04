import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Prepare scripts for square segmentations."
)
parser.add_argument("width", type=int, help="patch width.")
parser.add_argument("unit", choices=["pixel", "microns"], help="unit of measure.")
parser.add_argument("overlap", type=int, help="patch overlap.")
parser.add_argument(
    "-ir", "--intens_rat", default=0.1, type=float, help="intensity ratio."
)
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_rastered_{args.size}{args.unit}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_rastered_{args.size}{args.unit}/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 08:00:00
#SBATCH --mem=128G
#SBATCH -J rastered{args.size}_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/rastered{args.size}_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/rastered{args.size}_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa.sqsh"

mamba activate sopa
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/rastered_segmentation.py \
 {value} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Negative_Control_Rastered_{args.size} \
 {args.size} {args.overlap} {args.unit} {args.intens_rat}
""")
    f.close()
