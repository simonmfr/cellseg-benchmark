#!/bin/bash
# Submit cell type annotation of all methods of one sample; Negative_Control_Rastered_5 as separate job with more resources.
#   bash annotation_submit_sample.sh foxf2_s2_r1
set -euo pipefail

SAMPLE="${1:?usage: annotation_submit_sample.sh <sample>}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"

SAMPLE="$SAMPLE" sbatch "$DIR/annotation_one_sample_all_methods.sbatch"
if [[ -d "$ROOT/samples/$SAMPLE/results/Negative_Control_Rastered_5/sdata.zarr" ]]; then
  SAMPLE="$SAMPLE" SEG_METHOD='^Negative_Control_Rastered_5$' \
    sbatch -t 1-12:00:00 --mem=100G "$DIR/annotation_one_sample_one_method.sbatch"
fi
