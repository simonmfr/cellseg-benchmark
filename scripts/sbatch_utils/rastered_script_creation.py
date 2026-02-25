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
INPUT_PATH="{value["path"]}"

WIDTH="{args.width}"
OVERLAP="{args.overlap}"
UNIT="{args.unit}"
INTENS_RAT="{args.intens_rat}"

RESULT_DIR="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Negative_Control_Rastered_{args.width}"

CMD="python ~/gitrepos/cellseg-benchmark/scripts/segmentation/rastered_segmentation.py \\"${{INPUT_PATH}}\\" \\"${{RESULT_DIR}}\\" ${{WIDTH}} ${{OVERLAP}} ${{UNIT}} ${{INTENS_RAT}}"

write_log() {{
  local rc="$1"
  local end_iso="$2"
  local elapsed_s="$3"

  (
    flock -x 200
    if [ ! -f "${{RUN_LOG}}" ]; then
      printf "start_iso\tend_iso\telapsed_s\trc\tjobid\tjobname\tkey\tcp_version\tstaining\tconfidence\tinput_path\tresult_dir\thost\tnodelist\tsubmit_dir\tcmd\n" >> "${{RUN_LOG}}"
    fi
    # not a Cellpose/Proseg confidence-style run -> NA
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${{START_ISO}}" "${{end_iso}}" "${{elapsed_s}}" "${{rc}}" \
      "${{JOBID}}" "${{JOBNAME}}" "${{KEY}}" "NA" "NA" "NA" \
      "${{INPUT_PATH}}" "${{RESULT_DIR}}" "${{HOST}}" "${{NODELIST}}" "${{SUBMIT_DIR}}" "${{CMD}}" \
      >> "${{RUN_LOG}}"
  ) 200>>"${{LOCK_FILE}}"
}}

trap 'rc=$?; end_iso="$(date -Is)"; end_epoch="$(date +%s)"; elapsed_s=$((end_epoch-START_EPOCH)); write_log "$rc" "$end_iso" "$elapsed_s"' EXIT

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python ~/gitrepos/cellseg-benchmark/scripts/segmentation/rastered_segmentation.py \\
  "${{INPUT_PATH}}" \\
  "${{RESULT_DIR}}" \\
  "${{WIDTH}}" "${{OVERLAP}}" "${{UNIT}}" "${{INTENS_RAT}}"
""")
    f.close()
