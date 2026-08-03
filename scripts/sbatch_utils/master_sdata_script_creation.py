#!/usr/bin/env python
import pathlib
import shlex
import yaml

BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
YAML = BASE_PATH / "misc/sample_metadata.yaml"
OUT = BASE_PATH / "misc/sbatches/sbatch_master_sdata"
MANDATORY = {"cohort", "slide", "region", "organism", "run_date", "path"}

with open(YAML) as f:
    data = yaml.safe_load(f)

OUT.mkdir(parents=False, exist_ok=True)

for sample, meta in data.items():
    # extra obs = non-mandatory keys
    extras = []
    for k, v in meta.items():
        if k in MANDATORY or v is None or isinstance(v, (list, dict)):
            continue
        extras += ["--obs", f"{k}={v}"]

    argv = [
        sample,
        meta["path"],
        "z3",
        BASE_PATH,
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
    cli_args = " \\\n".join(shlex.quote(str(a)) for a in argv)

    sbatch = f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH -J master_sdata_{sample}
#SBATCH -o {BASE_PATH}/misc/logs/merged/%x.log
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -eu

mamba activate seg_postprocessing
python $HOME/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/master_sdata.py {cli_args}
"""
    (OUT / f"{sample}.sbatch").write_text(sbatch)
