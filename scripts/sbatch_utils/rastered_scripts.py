import argparse
from pathlib import Path

import yaml

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
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_rastered_{args.width}{args.unit}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_rastered_{args.width}{args.unit}/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 08:00:00
#SBATCH --mem=128G
#SBATCH -J rastered{args.width}_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/rastered{args.width}_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/rastered{args.width}_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/sopa_rastered_ABC.sqsh"

mamba activate rastered
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/segmentation/rastered_segmentation.py {value["path"]} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Negative_Control_Rastered_{args.width} \
 {args.width} {args.overlap} {args.unit} {args.intens_rat}
""")
    f.close()
