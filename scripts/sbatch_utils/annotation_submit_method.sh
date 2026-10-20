#!/bin/bash
# Submit cell type annotation of one method across all samples that still lack an annotation.
# Scans $ROOT/samples/*/results/<method>/sdata.zarr; a sample is "annotated" once
# cell_type_annotation/adata_obs_annotated.csv exists (last file written by cell_type_annotation.py;
# the directory itself is created at script start, so it is not a reliable marker).
# Negative_Control_Rastered_5 is submitted with more resources.
#
#   bash annotation_submit_method.sh Cellpose_1_nuclei_model
#   bash annotation_submit_method.sh Cellpose_1_nuclei_model --dry-run   # only list what would be queued
#   bash annotation_submit_method.sh Cellpose_1_nuclei_model --force     # also requeue annotated samples
set -euo pipefail

METHOD="${1:?usage: annotation_submit_method.sh <method> [--dry-run] [--force]}"
shift || true
DRY=0; FORCE=0
for a in "$@"; do
  case "$a" in
    --dry-run) DRY=1 ;;
    --force)   FORCE=1 ;;
    *) echo "unknown option: $a" >&2; exit 1 ;;
  esac
done

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
MARKER="cell_type_annotation/adata_obs_annotated.csv"

# exact-match regex for the one-method sbatch (escape regex metacharacters in the name)
ESCAPED="$(printf '%s' "$METHOD" | sed -e 's/[][\.^$*+?(){}|/]/\\&/g')"
SEG_REGEX="^${ESCAPED}$"

n_found=0; n_done=0; n_queued=0
shopt -s nullglob
for zarr in "$ROOT"/samples/*/results/"$METHOD"/sdata.zarr; do
  [[ -d "$zarr" ]] || continue
  mdir="$(dirname "$zarr")"
  SAMPLE="$(basename "$(dirname "$(dirname "$mdir")")")"
  n_found=$((n_found + 1))

  if (( ! FORCE )) && [[ -s "$mdir/$MARKER" ]]; then
    n_done=$((n_done + 1))
    continue
  fi

  if [[ "$METHOD" == "Negative_Control_Rastered_5" ]]; then
    OPTS=(-t 1-12:00:00 --mem=100G)
  else
    OPTS=()   # sbatch defaults: 06:00:00, 60G, 5 CPUs
  fi

  if (( DRY )); then
    echo "[dry-run] $SAMPLE | $METHOD ${OPTS[*]:-}"
  else
    SAMPLE="$SAMPLE" SEG_METHOD="$SEG_REGEX" \
      sbatch "${OPTS[@]}" "$DIR/annotation_one_sample_one_method.sbatch"
  fi
  n_queued=$((n_queued + 1))
done

echo "[$METHOD] sdata.zarr found: $n_found | already annotated: $n_done | queued: $n_queued$( ((DRY)) && echo ' (dry-run)')"
(( n_found > 0 )) || echo "[warn] no samples with $METHOD/sdata.zarr under $ROOT/samples" >&2
