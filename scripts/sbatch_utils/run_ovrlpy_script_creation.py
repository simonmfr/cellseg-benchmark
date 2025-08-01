import argparse
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Create scripts to run ovrlpy over one cohort."
)
parser.add_argument("cohort", help="Cohort name, e.g., 'foxf2'")
parser.add_argument(
    "--age",
    default=False,
    action="store_true",
    help="Age information is available. Otherwise, assume age: 6m.",
)
parser.add_argument(
    "--genotype",
    default=False,
    action="store_true",
    help="Genotype information is available. Otherwise, assume WT.",
)
args = parser.parse_args()

base_path = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
sbatch_path = f"{base_path}/misc/sbatches/sbatch_run_ovrlpy"
container_image = f"{base_path}/misc/cellseg_benchmark_2.sqsh"
log_path = f"{base_path}/misc/logs/merged"

Path(sbatch_path).mkdir(parents=False, exist_ok=True)


with open(f"{sbatch_path}/{args.cohort}_.sbatch", "w") as f:
    f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 05:00:00
#SBATCH --mem=100G
#SBATCH -J ovrlpy_{args.cohort}
#SBATCH -o {log_path}/%x.log
#SBATCH --container-image="{container_image}"
cd ~/gitrepos/cellseg-benchmark
git pull
mamba activate cellseg_benchmark
python scripts/run_ovrlpy.py {args.cohort}{" --age" if args.cohort == "aging" else ""}{" --genotype" if args.genotype else ""}
""")
