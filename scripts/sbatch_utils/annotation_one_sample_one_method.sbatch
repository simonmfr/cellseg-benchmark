#!/bin/bash
#
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 06:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=5
#SBATCH -J annot_single
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/merged/annot_single_%j.log
#SBATCH -D /home/ubuntu/gitrepos/cellseg-benchmark
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/annotation.sqsh"
#SBATCH --container-mounts=/dev:/dev

# Cell type annotation of one method of one sample. SEG_METHOD is a regex that must match exactly one method.
#
#   SAMPLE=foxf2_s2_r1 SEG_METHOD='^Cellpose_1_nuclei_model$' sbatch annotation_one_sample_one_method.sbatch
#   SAMPLE=foxf2_s2_r1 SEG_METHOD='^Negative_Control_Rastered_5$' sbatch -t 1-12:00:00 --mem=100G annotation_one_sample_one_method.sbatch

set -euo pipefail

SAMPLE="${SAMPLE:?set SAMPLE, e.g. SAMPLE=foxf2_s2_r1}"
SEG_METHOD="${SEG_METHOD:?set SEG_METHOD regex, e.g. SEG_METHOD='^Cellpose_1_nuclei_model$'}"

ROOT="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
RESULTS_DIR="$ROOT/samples/$SAMPLE/results"
[[ -d "$RESULTS_DIR" ]] || { echo "no results dir: $RESULTS_DIR"; exit 1; }

declare -a METHODS=()
shopt -s nocasematch
while IFS= read -r m; do
  [[ -d "$RESULTS_DIR/$m/sdata.zarr" ]] || continue
  [[ "$m" =~ ${SEG_METHOD} ]] && METHODS+=("$m")
done < <(find "$RESULTS_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
shopt -u nocasematch

NM=${#METHODS[@]}
if (( NM != 1 )); then
  echo "[error] $SAMPLE: SEG_METHOD='$SEG_METHOD' matched $NM methods with sdata.zarr (need exactly 1)."
  (( NM == 0 )) || printf '  - %s\n' "${METHODS[@]}"
  exit 1
fi

METHOD="${METHODS[0]}"
LOGFILE="$ROOT/misc/logs/merged/annot_${SAMPLE}_${METHOD}.log"
echo "running SAMPLE=$SAMPLE | METHOD=$METHOD -> $LOGFILE"
exec >"$LOGFILE" 2>&1

mamba activate annotation
python scripts/seg_postprocessing/cell_type_annotation.py "$SAMPLE" "$METHOD" "$ROOT"
