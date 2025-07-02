################################################################################
# 0. SCRIPT USAGE AND INPUT VALIDATION
################################################################################
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <SAMPLE_ID> /path/to/your/detected_transcripts.csv"
    echo "Example: $0 foxf2_s2_r1 /path/to/foxf2_s2_r1/detected_transcripts.csv"
    exit 1
fi

SAMPLE_ID_ARG="$1"
DETECTED_TRANSCRIPTS_INPUT_FILE_ARG="$2"

if [ ! -f "${DETECTED_TRANSCRIPTS_INPUT_FILE_ARG}" ]; then
    echo "ERROR: Input detected_transcripts.csv file not found: ${DETECTED_TRANSCRIPTS_INPUT_FILE_ARG}"
    exit 1
fi

################################################################################
# 1. ENVIRONMENT SETUP
################################################################################
echo "Starting script at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Hostname: $(hostname)"
echo "Sample ID: ${SAMPLE_ID_ARG}"
echo "Input detected_transcripts.csv: ${DETECTED_TRANSCRIPTS_INPUT_FILE_ARG}"

echo "Loading modules..."
module load slurm_setup

echo "Initializing Conda shell support..."
CONDA_BASE_DIR=$(conda info --base)
if [ -f "${CONDA_BASE_DIR}/etc/profile.d/conda.sh" ]; then
    source "${CONDA_BASE_DIR}/etc/profile.d/conda.sh"
else
    echo "ERROR: conda.sh not found in conda base directory: ${CONDA_BASE_DIR}"
    exit 1
fi

echo "Activating conda environment..."
conda activate ficture
# IMPORTANT: Ensure 'mpmath' package is installed in this environment for de_bulk:
# pip install mpmath

if [ "$CONDA_DEFAULT_ENV" != "ficture" ]; then
    echo "ERROR: Conda environment 'ficture' not activated correctly."
    echo "Current CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV"
    exit 1
fi
echo "Conda environment activated: $CONDA_DEFAULT_ENV"
echo "Python path: $(which python)"
echo "Ficture path: $(which ficture)"

################################################################################
# 2. PATH AND PARAMETER DEFINITIONS
################################################################################
echo "Defining paths and parameters..."

# --- Shared Resources Paths (Fixed) ---
LRZ_FICTURE_RESOURCES_DIR="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/Ficture/"
FORMAT_VIZGEN_SCRIPT="${LRZ_FICTURE_RESOURCES_DIR}/format_vizgen.py"
AGGREGATED_MODEL_PATH="${LRZ_FICTURE_RESOURCES_DIR}/aggregated_counts.tsv.gz"

# --- FICTURE Parameters ---
train_width=6 # Used for model_id and iden suffix
nFactor=21    # Used for model_id and iden suffix
thread=48     # Number of threads for ficture tools

# --- Define Output Structure ---
BASE_PROJECT_OUTPUT_DIR="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/"
MODEL_ID_PREFIX="nF${nFactor}.d_${train_width}" # e.g., nF21.d_6
OUTPUT_DIR_SUFFIX="-bulkRNAseq-noMBP" # Your specified suffix
FINAL_OUTPUT_DIR_NAME="${MODEL_ID_PREFIX}${OUTPUT_DIR_SUFFIX}" # e.g., nF21.d_6-bulkRNAseq-noMBP

# This is the root directory for all outputs of this specific run
RUN_SPECIFIC_OUTPUT_ROOT="${BASE_PROJECT_OUTPUT_DIR}/${SAMPLE_ID_ARG}/results/Ficture/output/${FINAL_OUTPUT_DIR_NAME}"
mkdir -p "${RUN_SPECIFIC_OUTPUT_ROOT}"
echo "Run-specific output root: ${RUN_SPECIFIC_OUTPUT_ROOT}"

# --- Path for NOMBP related outputs (Data preparation steps) ---
# This 'path' variable matches the usage in your original script for NOMBP files
path="${RUN_SPECIFIC_OUTPUT_ROOT}/NOMBP"
mkdir -p "${path}"
echo "NOMBP output path: ${path}"

# --- Path for analysis related outputs (Model fitting, DE, figures) ---
# The 'model_id' (e.g., nF21.d_6) will be used for file naming within this analysis path
# and also for subfolder naming if needed (as in original script's output_path_analysis_dir).
# For simplicity, 'analysis' itself can be the direct container for model_id named files/folders.
path_analysis_outputs="${RUN_SPECIFIC_OUTPUT_ROOT}/analysis"
mkdir -p "${path_analysis_outputs}"
echo "Analysis output path: ${path_analysis_outputs}"

# --- Identifier for filenames (using the user-provided SAMPLE_ID_ARG) ---
iden="${SAMPLE_ID_ARG}_width${train_width}_factor${nFactor}" # Changed to use SAMPLE_ID_ARG
echo "Run identifier (iden) for filenames: ${iden}"

################################################################################
# 3. FICTURE PIPELINE STEPS
################################################################################

# --- Step 1: Filtering detected_transcripts.csv ---
echo "Step 1: Filtering ${DETECTED_TRANSCRIPTS_INPUT_FILE_ARG}..."
filtered_transcripts_output_file="${path}/detected_transcripts_filtered.csv"
echo "Output of filtering: ${filtered_transcripts_output_file}"

grep -v "Mbp" "${DETECTED_TRANSCRIPTS_INPUT_FILE_ARG}" > "${filtered_transcripts_output_file}"
echo "Filtering done."

# --- Step 2: Formatting with format_vizgen.py ---
echo "Step 2: Formatting with format_vizgen.py..."
input_vizgen=${filtered_transcripts_output_file}
output_vizgen_matrix=${path}/filtered.matrix.${iden}.tsv
feature_vizgen=${path}/feature.clean.${iden}.tsv.gz
coor_vizgen=${path}/coordinate_minmax.tsv

python "${FORMAT_VIZGEN_SCRIPT}" \
    --input "${input_vizgen}" \
    --output "${output_vizgen_matrix}" \
    --feature "${feature_vizgen}" \
    --coor_minmax "${coor_vizgen}" \
    --precision 2 \
    --dummy_genes Blank
echo "format_vizgen.py done."

echo "Sorting and gzipping vizgen output..."
sort -S 4G -k1,1g "${output_vizgen_matrix}" | gzip -c > "${output_vizgen_matrix}.gz"
rm "${output_vizgen_matrix}"
echo "Sorting and gzipping vizgen output done."

# --- Step 3: FICTURE make_spatial_minibatch ---
echo "Step 3: FICTURE make_spatial_minibatch..."
mu_scale=1
key="Count"
major_axis="X"
batch_size=500
batch_buff=30
input_minibatch="${path}/filtered.matrix.${iden}.tsv.gz"
output_minibatch_base="${path}/batched.matrix.tsv"
output_minibatch_gz="${path}/batched.matrix.tsv.gz"

ficture make_spatial_minibatch --input "${input_minibatch}" --output "${output_minibatch_base}" --mu_scale ${mu_scale} --batch_size ${batch_size} --batch_buff ${batch_buff} --major_axis "${major_axis}"
echo "make_spatial_minibatch done."

echo "Sorting and gzipping minibatch output..."
sort -S 4G -k2,2n -k1,1g "${output_minibatch_base}" | gzip -c > "${output_minibatch_gz}"
rm "${output_minibatch_base}"
echo "Minibatch sorting and gzipping done."

# --- Step 4: FICTURE make_dge ---
echo "Step 4: FICTURE make_dge..."
# train_width is already defined
min_ct_per_unit=50
input_dge="${path}/filtered.matrix.${iden}.tsv.gz"
out_dge_base="${path}/hexagon.d_${train_width}.tsv"
out_dge_gz="${path}/hexagon.d_${train_width}.tsv.gz"

ficture make_dge --key "${key}" --count_header "${key}" --input "${input_dge}" --output "${out_dge_base}" --hex_width ${train_width} --n_move 2 --min_ct_per_unit ${min_ct_per_unit} --mu_scale ${mu_scale} --precision 2 --major_axis "${major_axis}"
echo "make_dge done."

echo "Sorting and gzipping DGE output..."
sort -S 4G -k1,1n "${out_dge_base}" | gzip -c > "${out_dge_gz}"
rm "${out_dge_base}"
echo "DGE sorting and gzipping done."

# --- Step 5: FICTURE Model Initialization/Fitting ---
echo "Step 5: FICTURE Model Initialization/Fitting..."
# nFactor, train_width, thread are already defined
model_id="${MODEL_ID_PREFIX}" # e.g., nF21.d_6 (This is for naming files and subdirs related to this model)
min_ct_per_feature=20

# Path for model-specific analysis outputs
output_path_model_analysis_dir="${path_analysis_outputs}/${model_id}" # e.g., .../analysis/nF21.d_6
mkdir -p "${output_path_model_analysis_dir}/figure/sub"

output_id=${model_id} # Consistent with original script's use of output_id for model files
output_fit_prefix="${output_path_model_analysis_dir}/${output_id}" # e.g., .../analysis/nF21.d_6/nF21.d_6
figure_path="${output_path_model_analysis_dir}/figure"          # e.g., .../analysis/nF21.d_6/figure

hexagon_input_fit="${path}/hexagon.d_${train_width}.tsv.gz"    # From NOMBP
feature_input_fit="${path}/feature.clean.${iden}.tsv.gz"       # From NOMBP

echo "Running ficture init_model_from_pseudobulk..."
ficture init_model_from_pseudobulk \
    --input "${hexagon_input_fit}" \
    --output "${output_fit_prefix}" \
    --feature "${feature_input_fit}" \
    --epoch 0 \
    --scale_model_rel -1 \
    --key "${key}" \
    --min_ct_per_feature ${min_ct_per_feature} \
    --thread ${thread} \
    --model "${AGGREGATED_MODEL_PATH}"
echo "init_model_from_pseudobulk done."

echo "Copying posterior count to model matrix..."
cp "${output_fit_prefix}.posterior.count.tsv.gz" "${output_fit_prefix}.model_matrix.tsv.gz"
echo "Copy done."

# --- Step 6: Visualization Preparation ---
echo "Step 6: Visualization Preparation..."
cmap_name="turbo"
input_choose_color="${output_fit_prefix}.fit_result.tsv.gz" # Output from init_model
cmap_path="${figure_path}/${output_id}.rgb.tsv"             # output_id is model_id here

echo "Running ficture choose_color..."
ficture choose_color --input "${input_choose_color}" --output "${figure_path}/${output_id}" --cmap_name "${cmap_name}"
echo "choose_color done."

input_plot_base="${output_fit_prefix}.fit_result.tsv.gz"
output_plot_base_prefix="${figure_path}/${output_id}.coarse"
fillr=$((train_width/2+1))

echo "Running ficture plot_base..."
ficture plot_base --input "${input_plot_base}" --output "${output_plot_base_prefix}" --fill_range ${fillr} --color_table "${cmap_path}" --plot_um_per_pixel 1 --plot_discretized
echo "plot_base done."

# --- Step 7: FICTURE Transform ---
echo "Step 7: FICTURE Transform..."
min_ct_per_unit_fit=20
fit_width=6 # Should match train_width for this model_id
anchor_res=4
fit_nmove=$((fit_width/anchor_res))
anchor_info="prj_${fit_width}.r_${anchor_res}"

# output_transform_prefix is based on output_id (model_id) and anchor_info, placed in model's analysis dir
output_transform_prefix="${output_path_model_analysis_dir}/${output_id}.${anchor_info}"
model_matrix_path="${output_fit_prefix}.model_matrix.tsv.gz"
pixel_input_transform="${path}/filtered.matrix.${iden}.tsv.gz" # From NOMBP

ficture transform --input "${pixel_input_transform}" --output_pref "${output_transform_prefix}" --model "${model_matrix_path}" --key "${key}" --major_axis "${major_axis}" --hex_width ${fit_width} --n_move ${fit_nmove} --min_ct_per_unit ${min_ct_per_unit_fit} --mu_scale ${mu_scale} --thread ${thread} --precision 2
echo "transform done."

# --- Step 8: Pixel level decoding & visualization (SLDA Decode) ---
echo "Step 8: Pixel level decoding (SLDA Decode)..."
radius=$(($anchor_res+1))
# prefix_decode uses output_id (model_id) and anchor_info
prefix_decode="${output_id}.${anchor_info}_${radius}"
input_slda_main="${path}/batched.matrix.tsv.gz" # From NOMBP
anchor_slda="${output_path_model_analysis_dir}/${output_id}.${anchor_info}.fit_result.tsv.gz" # Output from transform
output_slda_prefix="${output_path_model_analysis_dir}/${prefix_decode}" # Output slda results here
topk=3

ficture slda_decode --input "${input_slda_main}" --output "${output_slda_prefix}" --model "${model_matrix_path}" --anchor "${anchor_slda}" --anchor_in_um --neighbor_radius ${radius} --mu_scale ${mu_scale} --key "${key}" --precision 0.1 --lite_topk_output_pixel ${topk} --lite_topk_output_anchor ${topk} --thread ${thread}
echo "slda_decode done."

# --- Step 9: Formatting SLDA Decode Output ---
echo "Step 9: Formatting SLDA Decode Output..."
input_slda_formatted="${output_path_model_analysis_dir}/${prefix_decode}.pixel.tsv.gz"
output_slda_formatted_sorted="${output_path_model_analysis_dir}/${prefix_decode}.pixel.sorted.tsv.gz"

if [ ! -f "${input_slda_formatted}" ]; then
    echo "WARNING: Input for SLDA formatting not found: ${input_slda_formatted}"
else
    K=$( echo "$model_id" | sed 's/nF\([0-9]\{1,\}\)\..*/\1/' ) # model_id is nF<X>.d_<Y>

    if [ ! -f "${coor_vizgen}" ]; then # coor_vizgen is ${path}/coordinate_minmax.tsv (NOMBP)
        echo "ERROR: Coordinate file not found: ${coor_vizgen}"
        exit 1
    fi
    # Read coordinates from coor_vizgen
    # Ensure xmin, xmax, ymin, ymax are exported or use a different method to pass them to perl
    eval $(awk -F'\t' '{print $1 "=" $2}' "${coor_vizgen}")
    echo "Coordinates from file: xmin=${xmin}, xmax=${xmax}; ymin=${ymin}, ymax=${ymax}"


    if [ -z "${xmin}" ] || [ -z "${xmax}" ] || [ -z "${ymin}" ] || [ -z "${ymax}" ]; then
        echo "ERROR: Failed to load coordinates from ${coor_vizgen}"
        exit 1
    fi

    offsetx=${xmin}
    offsety=${ymin}
    rangex=$( echo "(${xmax} - ${xmin} + 0.5)/1+1" | bc )
    rangey=$( echo "(${ymax} - ${ymin} + 0.5)/1+1" | bc )
    bsize=2000
    scale=100
    header="##K=${K};TOPK=${topk}\n##BLOCK_SIZE=${bsize};BLOCK_AXIS=X;INDEX_AXIS=Y\n##OFFSET_X=${offsetx};OFFSET_Y=${offsety};SIZE_X=${rangex};SIZE_Y=${rangey};SCALE=${scale}\n#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3"

    BGZIP_CMD="$(which bgzip)"
    TABIX_CMD="$(which tabix)"

    (echo -e "${header}" && zcat "${input_slda_formatted}" | tail -n +2 | perl -slane -MMath::Round=nearest '$F[0]=int((nearest(0.1,$F[1])-$offx)/$bsize) * $bsize; $F[1]=int((nearest(0.1,$F[1])-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int((nearest(0.1,$F[2])-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- -bsize=${bsize} -scale=${scale} -offx=${offsetx} -offy=${offsety} | sort -S 4G -k1,1g -k3,3g ) | "${BGZIP_CMD}" -c > "${output_slda_formatted_sorted}"

    "${TABIX_CMD}" -f -s1 -b3 -e3 "${output_slda_formatted_sorted}"
    echo "SLDA output formatting and indexing done."
fi

# --- Step 10: Differential Expression (DE) ---
echo "Step 10: Differential Expression..."
max_pval_output=1e-3
min_fold_output=1.5
# input_de uses output_fit_prefix from model fitting step
input_de="${output_fit_prefix}.posterior.count.tsv.gz"
# output_de is based on prefix_decode from slda step, placed in model's analysis dir
output_de="${output_path_model_analysis_dir}/${prefix_decode}.bulk_chisq.tsv"

if [ -f "${input_de}" ]; then
    ficture de_bulk --input "${input_de}" --output "${output_de}" --min_ct_per_feature ${min_ct_per_feature} --max_pval_output ${max_pval_output} --min_fold_output ${min_fold_output} --thread ${thread}
    echo "DE analysis done."
else
    echo "WARNING: Input for DE analysis not found: ${input_de}"
fi

# --- Step 11: Reporting and Final Plotting ---
echo "Step 11: Reporting and Final Plotting..."
# output_report uses prefix_decode, placed in model's analysis dir
output_report="${output_path_model_analysis_dir}/${prefix_decode}.factor.info.html"

# Prepare posterior.count.tsv.gz for factor_report (it expects files based on --pref)
source_posterior_count_for_report="${output_fit_prefix}.posterior.count.tsv.gz" # Original model posterior
dest_posterior_count_for_report="${output_path_model_analysis_dir}/${prefix_decode}.posterior.count.tsv.gz" # Expected by factor_report

if [ -f "$source_posterior_count_for_report" ]; then
    echo "Copying $source_posterior_count_for_report to $dest_posterior_count_for_report for factor_report."
    cp "$source_posterior_count_for_report" "$dest_posterior_count_for_report"

    # cmap_path is already defined: ${figure_path}/${output_id}.rgb.tsv
    ficture factor_report --path "${output_path_model_analysis_dir}" --pref "${prefix_decode}" --color_table "${cmap_path}"
    echo "Factor report done."
else
    echo "WARNING: Source posterior count file for factor_report not found: $source_posterior_count_for_report"
fi

# input_plot_pixel uses output_slda_formatted_sorted from slda step
input_plot_pixel="${output_slda_formatted_sorted}"
# output_plot_pixel_full uses figure_path (model's analysis figure dir) and prefix_decode
output_plot_pixel_full="${figure_path}/${prefix_decode}.pixel.png"

if [ -f "${input_plot_pixel}" ]; then
    ficture plot_pixel_full --input "${input_plot_pixel}" --color_table "${cmap_path}" --output "${output_plot_pixel_full}" --plot_um_per_pixel 0.5 --full
    echo "Full pixel plot done."

    ficture plot_pixel_single --input "${input_plot_pixel}" --output "${output_plot_pixel_full}" --plot_um_per_pixel 0.5 --full --all
    echo "Single pixel plots (if any generated by --all) done."
else
    echo "WARNING: Input for pixel plotting not found: ${input_plot_pixel}"
fi

################################################################################
# End of Script
################################################################################
echo "Script finished at $(date)"
