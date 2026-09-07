#!/usr/bin/env python
"""Generate the two sbatch scripts for BANKSY-based brain-region mapping.

Step 1 builds the pre-QC raster reference, step 2 clusters it. Paths are
expanded here rather than inside the job, because the container overrides HOME.
"""

import argparse
import pathlib

BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

parser = argparse.ArgumentParser(
    description="sbatch scripts for BANKSY region mapping."
)
parser.add_argument("cohort", help="Cohort name, e.g. 'aging'")
parser.add_argument("--seg_method", default="Negative_Control_Rastered_25")
parser.add_argument("--repo", default="~/gitrepos/cellseg-benchmark")
parser.add_argument(
    "--image", default=str(BASE_PATH / "misc/enroot_images/benchmark_new.sqsh")
)
parser.add_argument("--env", default="seg_postprocessing", help="env for step 1")
parser.add_argument(
    "--banksy_env", default="seg_postprocessing", help="env with rpy2 and R Banksy"
)
args = parser.parse_args()

repo = pathlib.Path(args.repo).expanduser().resolve()
method_dir = BASE_PATH / "analysis" / args.cohort / args.seg_method
sbatch_dir = BASE_PATH / "misc" / "sbatches" / "sbatch_banksy"
sbatch_dir.mkdir(parents=True, exist_ok=True)

TEMPLATE = """#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem={mem}
#SBATCH -t {time}
#SBATCH -J {job}
#SBATCH -o {base}/misc/logs/merged/%x.log
#SBATCH --container-image="{image}"

set -eu
export PYTHONPATH={repo}
export MPLCONFIGDIR={base}/misc/tmp/mpl
mkdir -p $MPLCONFIGDIR
mamba activate {env}

{cmd}
"""

jobs = {
    f"{args.cohort}_regions": {
        "mem": "150G",
        "time": "04:00:00",
        "env": args.env,
        "cmd": f"python {repo}/scripts/seg_postprocessing/raster_adata_for_regions.py"
        f" {args.cohort} --seg_method {args.seg_method}",
    },
    f"{args.cohort}_banksy": {
        "mem": "150G",
        "time": "12:00:00",
        "env": args.banksy_env,
        "cmd": f"python {repo}/scripts/seg_postprocessing/banksy_clustering.py"
        f" {args.cohort} {method_dir}/adatas/adata_regions.h5ad.gz {method_dir}",
    },
}

for job, spec in jobs.items():
    path = sbatch_dir / f"{job}.sbatch"
    path.write_text(
        TEMPLATE.format(base=BASE_PATH, image=args.image, repo=repo, job=job, **spec)
    )
    print(path)

print(f"\nREG=$(sbatch --parsable {sbatch_dir}/{args.cohort}_regions.sbatch)")
print(f"sbatch --dependency=afterok:$REG {sbatch_dir}/{args.cohort}_banksy.sbatch")
