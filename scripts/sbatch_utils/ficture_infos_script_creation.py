import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(
    description="scripts for computing ficture statistics."
)
parser.add_argument(
    "--recompute", action="store_true", help="Consider genotype differentiation"
)
args = parser.parse_args()

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

#SBATCH -p lrz-hgx-a100-80x4
#SBATCH --exclude=lrz-hgx-a100-002
#SBATCH -t 09:00:00
#SBATCH --mem=900G
#SBATCH --gres=gpu:1
#SBATCH -J ficture_stats_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/ficture_stats_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/ficture_stats_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/cellseg_benchmark.sqsh"

cd ~/gitrepos/cellseg-benchmark
git pull
mamba activate cellseg_benchmark
python scripts/ficture_infos.py {key} {value} {"--recompute" if args.recompute else ""}
""")
    f.close()
