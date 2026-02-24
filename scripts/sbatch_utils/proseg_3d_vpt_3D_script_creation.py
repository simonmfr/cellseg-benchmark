import argparse
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(
    description="Prepare scripts for ProSeg with prior segmentation."
)
parser.add_argument("vpt_flavor", choices=["nuclei", "PolyT", "PolyT_nuclei"], help="vpt flavor.")
parser.add_argument("vpt_dim", choices=["2D", "3D"], help="vpt dimension.")
parser.add_argument("--voxel", default=1, type=int, help="number of z-layers.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Proseg_3D_vpt{args.vpt_dim}_{args.vpt_flavor}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Proseg_3D_vpt{args.vpt_dim}_{args.vpt_flavor}/{key}_{args.voxel}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 18:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=30
#SBATCH -J Proseg_3D_{key}_vpt{args.vpt_dim}_{args.vpt_flavor}_vxl_{args.voxel}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/Proseg_3D_{key}_vpt{args.vpt_dim}_{args.vpt_flavor}_vxl_{args.voxel}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/Proseg_3D_{key}_vpt{args.vpt_dim}_{args.vpt_flavor}_vxl_{args.voxel}.err
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
VPT_DIM="{args.vpt_dim}"
VPT_FLAVOR="{args.vpt_flavor}"
VOXEL="{args.voxel}"
INPUT_PATH="{value["path"]}"

RESULT_DIR="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Proseg_3D_vpt{args.vpt_dim}_{args.vpt_flavor}"

MODEL_NAME="vpt_${{VPT_DIM}}_DAPI_${{VPT_FLAVOR}}"

CMD="python ~/gitrepos/cellseg-benchmark/scripts/segmentation/proseg_3D_vpt_3D.py \\"${{INPUT_PATH}}\\" ${{KEY}} ${{MODEL_NAME}} --voxel-layers ${{VOXEL}} --output-cell-polygon-layers cell-polygons-layers.geojson.gz"

write_log() {{
  local rc="$1"
  local end_iso="$2"
  local elapsed_s="$3"

  (
    flock -x 200
    if [ ! -f "${{RUN_LOG}}" ]; then
      printf "start_iso\tend_iso\telapsed_s\trc\tjobid\tjobname\tkey\tcp_version\tstaining\tconfidence\tinput_path\tresult_dir\thost\tnodelist\tsubmit_dir\tcmd\n" >> "${{RUN_LOG}}"
    fi
    # not a Cellpose run -> cp_version/staining/confidence as NA; model info is in cmd/result_dir
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
python ~/gitrepos/cellseg-benchmark/scripts/segmentation/proseg_3D_vpt_3D.py \\
  "${{INPUT_PATH}}" \\
  "${{KEY}}" \\
  "${{MODEL_NAME}}" \\
  --voxel-layers "${{VOXEL}}" \\
  --output-cell-polygon-layers cell-polygons-layers.geojson.gz
""")
    f.close()
