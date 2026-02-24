from pathlib import Path
import yaml

YAML = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
SBATCH_DIR = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_visium"
OUT_DIR = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{k}/results/Negative_Control_Visium"

with open(YAML) as f:
    data = yaml.safe_load(f)

Path(SBATCH_DIR).mkdir(exist_ok=True)

for k, v in data.items():
    out = OUT_DIR.format(k=k)

    text = f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 10:00:00
#SBATCH --mem=128G
#SBATCH -J voronoi_{k}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/voronoi_{k}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/voronoi_{k}.err
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

KEY="{k}"
INPUT_PATH="{v["path"]}"
RESULT_DIR="{out}"

CMD="python ~/gitrepos/cellseg-benchmark/scripts/segmentation/voronoi_segmentation.py \\"${{INPUT_PATH}}\\" \\"${{RESULT_DIR}}\\""

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
python ~/gitrepos/cellseg-benchmark/scripts/segmentation/voronoi_segmentation.py \\
  "${{INPUT_PATH}}" \\
  "${{RESULT_DIR}}"
"""

    Path(SBATCH_DIR, f"{k}.sbatch").write_text(text)