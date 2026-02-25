import shlex
from pathlib import Path

import yaml

BASE = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
YAML = f"{BASE}/misc/sample_metadata.yaml"
OUT = f"{BASE}/misc/sbatches/sbatch_master_sdata"
MANDATORY = {"cohort", "slide", "region", "organism", "run_date", "path"}

with open(YAML) as f:
    data = yaml.safe_load(f)

Path(OUT).mkdir(parents=False, exist_ok=True)

for sample, meta in data.items():
    # extra obs = non-mandatory keys
    extras = []
    for k, v in meta.items():
        if k in MANDATORY or v is None or isinstance(v, (list, dict)):
            continue
        extras += ["--obs", f"{k}={v}"]

    argv = [
        "python",
        "~/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/master_sdata.py",
        sample,
        meta["path"],
        "z3",
        BASE,
        "--cohort",
        meta["cohort"],
        "--slide",
        str(meta["slide"]),
        "--region",
        str(meta["region"]),
        "--organism",
        meta["organism"],
        "--run_date",
        str(meta["run_date"]),
        *extras,
    ]
    cli = " \\\n".join(shlex.quote(str(a)) for a in argv)

    sbatch = f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH -J master_sdata_{sample}
#SBATCH -o {BASE}/misc/logs/merged/%x.log
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/benchmark.sqsh"

set -eu

mamba activate seg_postprocessing
"$CONDA_PREFIX/bin/time" -v \
  -o "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/master_sdata_{sample}_$(date +%Y%m%d_%H%M%S).time" \
{cli}
"""
    (Path(OUT) / f"{sample}.sbatch").write_text(sbatch)
