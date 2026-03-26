import argparse
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(description="scripts for Baysor segmentation.")
parser.add_argument("staining", help="Staining of prior cellpose segmentation.")
parser.add_argument("CP_version", help="Cellpose version of prior  segmentation.")
parser.add_argument("confidence", help="Confidence of prior cellpose segmentation.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Baysor_CP{args.CP_version}_{args.staining}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    if args.staining == "nuclei":
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Baysor_CP{args.CP_version}_{args.staining}/{key}_{args.confidence}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=75G
#SBATCH -J Baysor_{key}_CP1_{args.staining}_{args.confidence}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/Baysor_{key}_CP1_{args.staining}_{args.confidence}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/Baysor_{key}_CP1_{args.staining}_{args.confidence}.err
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
CP_VERSION="{args.CP_version}"
STAINING="{args.staining}"
CONFIDENCE="{args.confidence}"
INPUT_PATH="{value["path"]}"

RESULT_DIR="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Baysor_2D_Cellpose_1_{args.staining}_model_{args.confidence}"

CMD="python ~/gitrepos/cellseg-benchmark/scripts/segmentation/baysor.py \\"${{INPUT_PATH}}\\" Cellpose_1_${{STAINING}}_model ${{CONFIDENCE}} ${{KEY}}"

write_log() {{
  local rc="$1"
  local end_iso="$2"
  local elapsed_s="$3"

  (
    flock -x 200
    if [ ! -f "${{RUN_LOG}}" ]; then
      printf "start_iso\\tend_iso\\telapsed_s\\trc\\tjobid\\tjobname\\tkey\\tcp_version\\tstaining\\tconfidence\\tinput_path\\tresult_dir\\thost\\tnodelist\\tsubmit_dir\\tcmd\\n" >> "${{RUN_LOG}}"
    fi
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
      "${{START_ISO}}" "${{end_iso}}" "${{elapsed_s}}" "${{rc}}" \\
      "${{JOBID}}" "${{JOBNAME}}" "${{KEY}}" "${{CP_VERSION}}" "${{STAINING}}" "${{CONFIDENCE}}" \\
      "${{INPUT_PATH}}" "${{RESULT_DIR}}" "${{HOST}}" "${{NODELIST}}" "${{SUBMIT_DIR}}" "${{CMD}}" \\
      >> "${{RUN_LOG}}"
  ) 200>>"${{LOCK_FILE}}"
}}

trap 'rc=$?; end_iso="$(date -Is)"; end_epoch="$(date +%s)"; elapsed_s=$((end_epoch-START_EPOCH)); write_log "$rc" "$end_iso" "$elapsed_s"' EXIT

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python ~/gitrepos/cellseg-benchmark/scripts/segmentation/baysor.py \
  "${{INPUT_PATH}}" \
  "Cellpose_1_${{STAINING}}"_model \
  "${{CONFIDENCE}}" \
  "${{KEY}}"
""")
        f.close()
    else:
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Baysor_CP{args.CP_version}_{args.staining}/{key}_{args.confidence}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=75G
#SBATCH -J Baysor_{key}_CP{args.CP_version}_{args.staining}_{args.confidence}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/Baysor_{key}_CP{args.CP_version}_{args.staining}_{args.confidence}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/Baysor_{key}_CP{args.CP_version}_{args.staining}_{args.confidence}.err
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
CP_VERSION="{args.CP_version}"
STAINING="{args.staining}"
CONFIDENCE="{args.confidence}"
INPUT_PATH="{value["path"]}"

RESULT_DIR="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Baysor_2D_Cellpose_{args.CP_version}_DAPI_{args.staining}_{args.confidence}"

CMD="python ~/gitrepos/cellseg-benchmark/scripts/segmentation/baysor.py \\"${{INPUT_PATH}}\\" Cellpose_${{CP_VERSION}}_DAPI_${{STAINING}} ${{CONFIDENCE}} ${{KEY}}"

write_log() {{
  local rc="$1"
  local end_iso="$2"
  local elapsed_s="$3"

  (
    flock -x 200
    if [ ! -f "${{RUN_LOG}}" ]; then
      printf "start_iso\\tend_iso\\telapsed_s\\trc\\tjobid\\tjobname\\tkey\\tcp_version\\tstaining\\tconfidence\\tinput_path\\tresult_dir\\thost\\tnodelist\\tsubmit_dir\\tcmd\\n" >> "${{RUN_LOG}}"
    fi
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
      "${{START_ISO}}" "${{end_iso}}" "${{elapsed_s}}" "${{rc}}" \\
      "${{JOBID}}" "${{JOBNAME}}" "${{KEY}}" "${{CP_VERSION}}" "${{STAINING}}" "${{CONFIDENCE}}" \\
      "${{INPUT_PATH}}" "${{RESULT_DIR}}" "${{HOST}}" "${{NODELIST}}" "${{SUBMIT_DIR}}" "${{CMD}}" \\
      >> "${{RUN_LOG}}"
  ) 200>>"${{LOCK_FILE}}"
}}

trap 'rc=$?; end_iso="$(date -Is)"; end_epoch="$(date +%s)"; elapsed_s=$((end_epoch-START_EPOCH)); write_log "$rc" "$end_iso" "$elapsed_s"' EXIT

mamba activate segmentation

mkdir -p "${{RESULT_DIR}}"
python ~/gitrepos/cellseg-benchmark/scripts/segmentation/baysor.py \
  "${{INPUT_PATH}}" \
  "Cellpose_${{CP_VERSION}}_DAPI_${{STAINING}}" \
  "${{CONFIDENCE}}" \
  "${{KEY}}"
""")
        f.close()
