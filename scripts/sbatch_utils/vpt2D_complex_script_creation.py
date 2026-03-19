import argparse
from os import listdir
from os.path import join
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(description="Prepare scripts for vpt pipeline.")
parser.add_argument("staining1", help="Name of cell-boundary staining.")
parser.add_argument("staining2", help="Name of nucleus staining.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

experiment_json_path = str(Path(__file__).parents[2] / "configs" / f"vpt_2D_{args.staining1}_{args.staining2}.json")

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_2D_complex"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    res_path = f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/vpt_2D_DAPI_{args.staining1}_{args.staining2}"
    for dire in listdir(value["path"]):
        if dire.endswith(".vzg"):
            vzg_path = join(value["path"], dire)
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_2D_complex/{key}_{args.staining1}_{args.staining2}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
#SBATCH -p lrz-hgx-a100-80x4
#SBATCH -t 1-12:00:00
#SBATCH --mem=600G
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=40
#SBATCH -J vpt2D_{key}_{args.staining1}_{args.staining2}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/vpt2D_{key}_{args.staining1}_{args.staining2}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/vpt2D_{key}_{args.staining1}_{args.staining2}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/vpt.sqsh"

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
STAINING1="{args.staining1}"
STAINING2="{args.staining2}"

RES_PATH="{res_path}"
EXPERIMENT_JSON="{experiment_json_path}"
INPUT_IMAGES="{join(value["path"], "images")}"
INPUT_TRANSFORM="{join(value["path"], "images/micron_to_mosaic_pixel_transform.csv")}"
INPUT_TRANSCRIPTS="{join(value["path"], "detected_transcripts.csv")}"
VZG_PATH="{vzg_path}"

BOUNDARIES="{join(res_path, "analysis_outputs/cellpose2_micron_space.parquet")}"
OUT_ANALYSIS="{join(res_path, "analysis_outputs")}"
OUT_CBG="{join(res_path, "analysis_outputs/cell_by_gene.csv")}"
OUT_META="{join(res_path, "analysis_outputs/cell_metadata.csv")}"
OUT_VZG="{join(res_path, "visualize.vzg")}"
TMP_PATH="{join(res_path, "tmp")}"

CMD="vpt run-segmentation (json=${{EXPERIMENT_JSON}}; images=${{INPUT_IMAGES}}; out=${{OUT_ANALYSIS}}) + partition-transcripts + derive-entity-metadata + update-vzg"

write_log() {{
  local rc="$1"
  local end_iso="$2"
  local elapsed_s="$3"

  (
    flock -x 200
    if [ ! -f "${{RUN_LOG}}" ]; then
      printf "start_iso\tend_iso\telapsed_s\trc\tjobid\tjobname\tkey\tcp_version\tstaining\tconfidence\tinput_path\tresult_dir\thost\tnodelist\tsubmit_dir\tcmd\n" >> "${{RUN_LOG}}"
    fi
    # not a Cellpose/Proseg confidence-style run -> NA; keep stainings in 'staining' field
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${{START_ISO}}" "${{end_iso}}" "${{elapsed_s}}" "${{rc}}" \
      "${{JOBID}}" "${{JOBNAME}}" "${{KEY}}" "NA" "${{STAINING1}}+${{STAINING2}}" "NA" \
      "${{INPUT_IMAGES}}" "${{RES_PATH}}" "${{HOST}}" "${{NODELIST}}" "${{SUBMIT_DIR}}" "${{CMD}}" \
      >> "${{RUN_LOG}}"
  ) 200>>"${{LOCK_FILE}}"
}}

trap 'rc=$?; end_iso="$(date -Is)"; end_epoch="$(date +%s)"; elapsed_s=$((end_epoch-START_EPOCH)); write_log "$rc" "$end_iso" "$elapsed_s"' EXIT

mamba activate vpt

mkdir -p "${{RES_PATH}}"

vpt --verbose --processes 40 run-segmentation \\
  --segmentation-algorithm "${{EXPERIMENT_JSON}}" \\
  --input-images "${{INPUT_IMAGES}}" \\
  --input-micron-to-mosaic "${{INPUT_TRANSFORM}}" \\
  --output-path "${{OUT_ANALYSIS}}" \\
  --tile-size 2400 \\
  --tile-overlap 200

vpt --verbose partition-transcripts \\
  --input-boundaries "${{BOUNDARIES}}" \\
  --input-transcripts "${{INPUT_TRANSCRIPTS}}" \\
  --output-entity-by-gene "${{OUT_CBG}}"

vpt --verbose derive-entity-metadata \\
  --input-boundaries "${{BOUNDARIES}}" \\
  --output-metadata "${{OUT_META}}"

vpt --verbose --processes 10 update-vzg \\
  --input-vzg "${{VZG_PATH}}" \\
  --input-boundaries "${{BOUNDARIES}}" \\
  --input-entity-by-gene "${{OUT_CBG}}" \\
  --output-vzg "${{OUT_VZG}}" \\
  --input-metadata "${{OUT_META}}" \\
  --temp-path "${{TMP_PATH}}"
""")
    f.close()
