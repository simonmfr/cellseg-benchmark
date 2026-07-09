#!/usr/bin/env bash
set -euo pipefail
export PATH=~/perl5/bin:$PATH
export PERL5LIB=~/perl5/lib/perl5:${PERL5LIB:-}

# --- Args ---
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <SAMPLE_ID> /path/to/detected_transcripts.csv"
    echo "Example: $0 foxf2_s2_r1 /path/to/foxf2_s2_r1/detected_transcripts.csv"
    exit 1
fi
SAMPLE_ID="$1"
DETECTED_TRANSCRIPTS="$2"
[ -f "$DETECTED_TRANSCRIPTS" ] || { echo "ERROR: input not found: $DETECTED_TRANSCRIPTS"; exit 1; }

# --- Helper: sort a file, gzip to <file>.gz, drop the original ---
sort_gzip() {  # usage: sort_gzip <file> <sort-args...>
    local f="$1"; shift
    sort -S 4G "$@" "$f" | gzip -c > "$f.gz"
    rm "$f"
}

# --- Env ---
echo "Start $(date) | job ${SLURM_JOB_ID:-NA} | host $(hostname) | sample $SAMPLE_ID"
#module load slurm_setup
#source "$(conda info --base)/etc/profile.d/conda.sh"
#conda activate ficture          # env must contain: mpmath (needed by de_bulk)
[ "$CONDA_DEFAULT_ENV" = "ficture" ] || { echo "ERROR: env 'ficture' not active"; exit 1; }
echo "python: $(which python) | ficture: $(which ficture)"

# --- Paths & parameters ---
# format_vizgen.py is loaded from the git repo (scripts/ficture/ next to this script).
# Override with FORMAT_VIZGEN_SCRIPT=/path/to/format_vizgen.py if needed.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FORMAT_VIZGEN_SCRIPT="${FORMAT_VIZGEN_SCRIPT:-${SCRIPT_DIR}/format_vizgen.py}"
[ -f "$FORMAT_VIZGEN_SCRIPT" ] || { echo "ERROR: format_vizgen.py not found: $FORMAT_VIZGEN_SCRIPT"; exit 1; }

RESOURCES_DIR="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc"
AGGREGATED_COUNTS="${RESOURCES_DIR}/scRNAseq_ref_ABCAtlas_Yao2023Nature/pseudobulk_celltype_for_ficture.norm.tsv.gz"

train_width=6
thread=32
mu_scale=1
key="Count"
major_axis="X"
min_ct_per_feature=20

BASE_OUT="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples"
RUN_ROOT="${BASE_OUT}/${SAMPLE_ID}/results/Ficture/output"
figure_path="${RUN_ROOT}/figure"
path="${RUN_ROOT}/tmp"                             # intermediates, deleted after a successful run
mkdir -p "$path" "$figure_path"

iden="${SAMPLE_ID}_width${train_width}"
fit_prefix="${RUN_ROOT}/model"                     # base path for fitted-model files
matrix="${path}/filtered.matrix.${iden}.tsv"
feature="${path}/feature.clean.${iden}.tsv.gz"
coor="${path}/coordinate_minmax.tsv"

# --- Step 1: filter out Mbp transcripts ---
echo "Step 1: filter Mbp"
filtered_csv="${path}/detected_transcripts_filtered.csv"
grep -v "Mbp" "$DETECTED_TRANSCRIPTS" > "$filtered_csv"

# --- Step 2: format with format_vizgen.py ---
echo "Step 2: format_vizgen.py"
python "$FORMAT_VIZGEN_SCRIPT" --input "$filtered_csv" --output "$matrix" \
    --feature "$feature" --coor_minmax "$coor" --precision 2 --dummy_genes Blank
sort_gzip "$matrix" -k1,1g
matrix_gz="${matrix}.gz"

# --- Step 3: make_spatial_minibatch ---
echo "Step 3: make_spatial_minibatch"
batched="${path}/batched.matrix.tsv"
ficture make_spatial_minibatch --input "$matrix_gz" --output "$batched" \
    --mu_scale $mu_scale --batch_size 500 --batch_buff 30 --major_axis "$major_axis"
sort_gzip "$batched" -k2,2n -k1,1g
batched_gz="${batched}.gz"

# --- Step 4: make_dge ---
echo "Step 4: make_dge"
hexagon="${path}/hexagon.d_${train_width}.tsv"
ficture make_dge --key "$key" --count_header "$key" --input "$matrix_gz" --output "$hexagon" \
    --hex_width $train_width --n_move 2 --min_ct_per_unit 50 --mu_scale $mu_scale \
    --precision 2 --major_axis "$major_axis"
sort_gzip "$hexagon" -k1,1n
hexagon_gz="${hexagon}.gz"

# --- Step 5: init model from pseudobulk ---
echo "Step 5: init_model_from_pseudobulk"
ficture init_model_from_pseudobulk --input "$hexagon_gz" --output "$fit_prefix" \
    --feature "$feature" --epoch 0 --scale_model_rel -1 --key "$key" \
    --min_ct_per_feature $min_ct_per_feature --thread $thread --model "$AGGREGATED_COUNTS"
cp "${fit_prefix}.posterior.count.tsv.gz" "${fit_prefix}.model_matrix.tsv.gz"
model_matrix="${fit_prefix}.model_matrix.tsv.gz"

# --- Step 6: load color + plot overview ---
echo "Step 6: choose color + plot overview"
GLOBAL_CMAP="${RESOURCES_DIR}/ficture_model.rgb.tsv"
cmap_path="${figure_path}/ficture_model.rgb.tsv"
if [ ! -f "$GLOBAL_CMAP" ]; then
  ficture choose_color --input "${fit_prefix}.fit_result.tsv.gz" \
    --output "${RESOURCES_DIR}/model.tmp.$$" --cmap_name turbo
  mv -n "${RESOURCES_DIR}/model.tmp.$$.rgb.tsv" "$GLOBAL_CMAP"
fi
cp "$GLOBAL_CMAP" "$cmap_path"
ficture plot_base --input "${fit_prefix}.fit_result.tsv.gz" \
    --output "${figure_path}/model.coarse" --fill_range $((train_width/2+1)) \
    --color_table "$cmap_path" --plot_um_per_pixel 1 --plot_discretized

# --- Step 7: transform ---
echo "Step 7: transform"
anchor_res=4
transform_prefix="${RUN_ROOT}/transform"
ficture transform --input "$matrix_gz" --output_pref "$transform_prefix" --model "$model_matrix" \
    --key "$key" --major_axis "$major_axis" --hex_width $train_width --n_move $((train_width/anchor_res)) \
    --min_ct_per_unit 20 --mu_scale $mu_scale --thread $thread --precision 2

# --- Step 8: pixel-level decoding (slda_decode) ---
echo "Step 8: slda_decode"
radius=$((anchor_res+1))
topk=3
prefix_decode="decode"
ficture slda_decode --input "$batched_gz" --output "${RUN_ROOT}/${prefix_decode}" \
    --model "$model_matrix" --anchor "${transform_prefix}.fit_result.tsv.gz" --anchor_in_um \
    --neighbor_radius $radius --mu_scale $mu_scale --key "$key" --precision 0.1 \
    --lite_topk_output_pixel $topk --lite_topk_output_anchor $topk --thread $thread

# --- Step 9: reformat + index slda output ---
echo "Step 9: format SLDA output"
pixel_in="${RUN_ROOT}/${prefix_decode}.pixel.tsv.gz"
pixel_sorted="${RUN_ROOT}/decode.pixel.sorted.tsv.gz"
if [ -f "$pixel_in" ]; then
    K=$(zcat "$model_matrix" | awk -F'\t' 'NR==1{print NF-1}')   # factor count = model matrix columns minus the gene column
    eval "$(awk -F'\t' '{print $1 "=" $2}' "$coor")"
    : "${xmin:?missing coords}" "${xmax:?missing coords}" "${ymin:?missing coords}" "${ymax:?missing coords}"
    echo "Coords: xmin=$xmin xmax=$xmax ymin=$ymin ymax=$ymax"
    offsetx=$xmin; offsety=$ymin
    rangex=$(echo "(${xmax} - ${xmin} + 0.5)/1+1" | bc)
    rangey=$(echo "(${ymax} - ${ymin} + 0.5)/1+1" | bc)
    bsize=2000; scale=100
    header="##K=${K};TOPK=${topk}\n##BLOCK_SIZE=${bsize};BLOCK_AXIS=X;INDEX_AXIS=Y\n##OFFSET_X=${offsetx};OFFSET_Y=${offsety};SIZE_X=${rangex};SIZE_Y=${rangey};SCALE=${scale}\n#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3"
    (echo -e "$header" && zcat "$pixel_in" | tail -n +2 | perl -MGetopt::Long -MMath::Round -lane 'BEGIN {GetOptions("bsize=f" => \$bsize, "scale=f" => \$scale, "offx=f"  => \$offx, "offy=f"  => \$offy,);}$F[0] = int((nearest(0.1,$F[1])-$offx)/$bsize) * $bsize; $F[1] = int((nearest(0.1,$F[1])-$offx)*$scale); $F[1] = ($F[1]>=0)?$F[1]:0; $F[2] = int((nearest(0.1,$F[2])-$offy)*$scale); $F[2] = ($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- --bsize=${bsize} --scale=${scale} --offx=${offsetx} --offy=${offsety} | sort -S 4G -k1,1g -k3,3g) | bgzip -c > "$pixel_sorted"
    tabix -f -s1 -b3 -e3 "$pixel_sorted"
else
    echo "WARNING: SLDA pixel input not found: $pixel_in"
fi

# --- Step 10: differential expression (de_bulk) ---
echo "Step 10: de_bulk"
de_in="${fit_prefix}.posterior.count.tsv.gz"
if [ -f "$de_in" ]; then
    ficture de_bulk --input "$de_in" --output "${RUN_ROOT}/${prefix_decode}.bulk_chisq.tsv" \
        --min_ct_per_feature $min_ct_per_feature --max_pval_output 1e-3 --min_fold_output 1.5 --thread $thread
else
    echo "WARNING: DE input not found: $de_in"
fi

# --- Step 11: factor report + pixel plots ---
echo "Step 11: factor_report + pixel plots"
if [ -f "${fit_prefix}.posterior.count.tsv.gz" ]; then
    cp "${fit_prefix}.posterior.count.tsv.gz" "${RUN_ROOT}/${prefix_decode}.posterior.count.tsv.gz"
    ficture factor_report --path "$RUN_ROOT" --pref "$prefix_decode" --color_table "$cmap_path"
else
    echo "WARNING: posterior count for factor_report not found"
fi

pixel_png="${figure_path}/${prefix_decode}.pixel.png"
if [ -f "$pixel_sorted" ]; then
    ficture plot_pixel_full --input "$pixel_sorted" --color_table "$cmap_path" \
        --output "$pixel_png" --plot_um_per_pixel 0.5 --full
    ficture plot_pixel_single --input "$pixel_sorted" --output "$pixel_png" \
        --plot_um_per_pixel 0.5 --full --all
else
    echo "WARNING: pixel plot input not found: $pixel_sorted"
fi

echo "Cleanup: removing $path and redundant files"
# redundant with model.posterior.count.tsv.gz
rm -f "${fit_prefix}.model_matrix.tsv.gz" "${RUN_ROOT}/${prefix_decode}.posterior.count.tsv.gz"
# superseded by the sorted + tabix-indexed file
[ -s "$pixel_sorted" ] && [ -s "${pixel_sorted}.tbi" ] && rm -f "$pixel_in"
rm -rf "$path"

echo "Done $(date)"