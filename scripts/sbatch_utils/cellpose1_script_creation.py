import argparse
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(description="scripts for Cellpose 1 segmentation.")
parser.add_argument("staining", help="Staining of prior cellpose segmentation.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Cellpose_1"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Cellpose_1/{key}_{args.staining}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH -J CP1_{key}_{args.staining}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/CP1_{key}_{args.staining}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/CP1_{key}_{args.staining}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/benchmark.sqsh"

set -euo pipefail

# ---------- central run log (shared across all scripts) ----------
RUN_LOG="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/job_runs.tsv"
LOCK_FILE="${{RUN_LOG}}.lock"
mkdir -p "$(dirname "${{RUN_LOG}}")"

JOBID="${{SLURM_JOB_ID:-NA}}"
JOBNAME="${{SLURM_JOB_NAME:-NA}}"
NODELIST="${{SLURM_JOB_NODELIST:-NA}}"
SUBMIT_DIR="${{SLURM_SUBMIT_DIR:-$PWD}}"
HOST="$(hostname -f 2>/dev/null || hostname)"

START_ISO="$(date -Is)"
START_EPOCH="$(date +%s)"

KEY="{key}"
STAINING="{args.staining}"
INPUT_PATH="{value["path"]}"

RESULT_DIR="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Cellpose_1_DAPI_{args.staining}"

CMD="python ~/gitrepos/cellseg-benchmark/scripts/segmentation/cellpose_1.py \\"${{INPUT_PATH}}\\" \\"${{RESULT_DIR}}\\" ${{STAINING}}"

write_log() {{
  local rc="$1"
  local end_iso="$2"
  local elapsed_s="$3"

  (
    flock -x 200
    if [ ! -f "${{RUN_LOG}}" ]; then
      printf "start_iso\\tend_iso\\telapsed_s\\trc\\tjobid\\tjobname\\tkey\\tcp_version\\tstaining\\tconfidence\\tinput_path\\tresult_dir\\thost\\tnodelist\\tsubmit_dir\\tcmd\\n" >> "${{RUN_LOG}}"
    fi
    # cp_version/confidence not applicable here -> NA
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
      "${{START_ISO}}" "${{end_iso}}" "${{elapsed_s}}" "${{rc}}" \\
      "${{JOBID}}" "${{JOBNAME}}" "${{KEY}}" "1" "${{STAINING}}" "NA" \\
      "${{INPUT_PATH}}" "${{RESULT_DIR}}" "${{HOST}}" "${{NODELIST}}" "${{SUBMIT_DIR}}" "${{CMD}}" \\
      >> "${{RUN_LOG}}"
  ) 200>>"${{LOCK_FILE}}"
}}

trap 'rc=$?; end_iso="$(date -Is)"; end_epoch="$(date +%s)"; elapsed_s=$((end_epoch-START_EPOCH)); write_log "$rc" "$end_iso" "$elapsed_s"' EXIT

mamba activate segmentation
mkdir -p "${{RESULT_DIR}}"
python ~/gitrepos/cellseg-benchmark/scripts/segmentation/cellpose_1.py \\
  "${{INPUT_PATH}}" \\
  "${{RESULT_DIR}}" \\
  "${{STAINING}}"
""")
    f.close()
