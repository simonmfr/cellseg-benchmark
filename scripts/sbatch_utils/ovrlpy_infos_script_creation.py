import argparse
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(description="scripts for computing ovrlpy statistics.")
parser.add_argument(
    "--recompute", action="store_true", help="Consider genotype differentiation"
)
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_ovrlpy_stats"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_ovrlpy_stats/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH -q cpu
#SBATCH -t 01:30:00
#SBATCH --mem=50G
#SBATCH -J ovrlpy_stats_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/ovrlpy_stats_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/ovrlpy_stats_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/cellseg_benchmark_2.sqsh"

cd ~/gitrepos/cellseg-benchmark
git pull
mamba activate cellseg_benchmark
python scripts/ovrlpy_infos.py {key} {value["path"]} {"--recompute" if args.recompute else ""}
""")
    f.close()
