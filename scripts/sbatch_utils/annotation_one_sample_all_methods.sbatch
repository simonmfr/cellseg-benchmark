#!/bin/bash
#
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 05:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=5
#SBATCH -J annot_sample
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/merged/annot_sample_%A_%a.log
#SBATCH --array=0-44%5
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/benchmark_new2.sqsh"

# Cell type annotation of all methods (sdata.zarr) of one sample, one array task per method.
# Negative_Control_Rastered_5 is skipped unless SEG_METHOD is set; annotation_submit_sample.sh
# submits it separately with more resources.
#
#   SAMPLE=foxf2_s2_r1 sbatch annotation_one_sample_all_methods.sbatch
#   SAMPLE=foxf2_s2_r1 SEG_METHOD='^(cellpose|baysor)' sbatch annotation_one_sample_all_methods.sbatch

set -euo pipefail

SAMPLE="${SAMPLE:?set SAMPLE, e.g. SAMPLE=foxf2_s2_r1}"

ROOT="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
RESULTS_DIR="$ROOT/samples/$SAMPLE/results"
[[ -d "$RESULTS_DIR" ]] || { echo "no results dir: $RESULTS_DIR"; exit 1; }

declare -a METHODS=()
shopt -s nocasematch
while IFS= read -r m; do
  [[ -d "$RESULTS_DIR/$m/sdata.zarr" ]] || continue
  if [[ -z "${SEG_METHOD:-}" ]]; then
    [[ "$m" == Negative_Control_Rastered_5 ]] || METHODS+=("$m")
  elif [[ "$m" =~ ${SEG_METHOD} ]]; then
    METHODS+=("$m")
  fi
done < <(find "$RESULTS_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
shopt -u nocasematch

NM=${#METHODS[@]}
if (( NM == 0 )); then
  echo "[skip] $SAMPLE: no methods with sdata.zarr${SEG_METHOD:+ matching SEG_METHOD='$SEG_METHOD'}"
  exit 0
fi
if (( NM > 45 )); then
  echo "[error] $SAMPLE has $NM methods (>45); narrow SEG_METHOD or raise --array (max 50 on lrz-cpu)."
  exit 1
fi
(( SLURM_ARRAY_TASK_ID < NM )) || exit 0

METHOD="${METHODS[$SLURM_ARRAY_TASK_ID]}"
LOGFILE="$ROOT/misc/logs/merged/annot_${SAMPLE}_${METHOD}.log"
echo "running SAMPLE=$SAMPLE | METHOD=$METHOD -> $LOGFILE"
exec >"$LOGFILE" 2>&1

mamba activate seg_postprocessing
python $HOME/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/cell_type_annotation.py "$SAMPLE" "$METHOD" "$ROOT"
