#!/usr/bin/env bash
#SBATCH --job-name=variant_effect_prediction
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:20:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ChromBPNet variant-effect prediction pipeline
#
# This script:
#   1. Scores variants with each model fold.
#   2. Averages variant scores across folds.
#   3. Annotates averaged variants with peaks, nearest genes, and FiNeMo hits.
#   4. Optionally computes per-fold variant SHAP values.
#
# No data are included with this script. Edit only the CONFIGURATION section
# below to match your environment, file layout, condition names, and models.
#
# Run with:
#   sbatch chrombpnet_variant_effect_pipeline.sh
#
# Optional environment overrides:
#   RUN_SHAP=true sbatch chrombpnet_variant_effect_pipeline.sh

set -Eeuo pipefail
IFS=$'\n\t'
trap 'echo "ERROR: command failed at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

################################################################################
# CONFIGURATION — EDIT THIS SECTION
################################################################################

# -----------------------------
# Cluster environment
# -----------------------------
# Set LOAD_MODULES=false if your environment is already configured.
LOAD_MODULES=true
MODULES=(
  "miniforge"
  "cuda/12.4.1"
)

# Conda environment containing ChromBPNet and its dependencies.
# Leave CONDA_INIT_SCRIPT empty to derive it from `conda info --base`.
CONDA_ENV="/path/to/chrombpnet/conda/environment"
CONDA_INIT_SCRIPT=""
PYTHON_BIN="python"

# -----------------------------
# Analysis labels and folds
# -----------------------------
CELL_TYPE="cell_type"
CONDITIONS=("condition_1" "condition_2")
FOLDS=(0 1 2 3 4)

# -----------------------------
# Project directories
# -----------------------------
PROJECT_DIR="/path/to/project"
MODEL_ROOT="${PROJECT_DIR}/trained_models"
OUTPUT_DIR="${PROJECT_DIR}/variant_effect_prediction"
ANNOTATION_INPUT_DIR="${PROJECT_DIR}/annotation_inputs"
HITS_SOURCE_DIR="${PROJECT_DIR}/modisco_reports"

# Template tokens supported below:
#   {cell_type}  {condition}  {fold}
MODEL_PATH_TEMPLATE="${MODEL_ROOT}/{cell_type}_{condition}_{fold}/models/{cell_type}_{condition}_{fold}_chrombpnet_nobias.h5"
HIT_FILE_GLOB_TEMPLATE="${HITS_SOURCE_DIR}/{cell_type}_{condition}_*/finemo_footprint/*hits_unique.tsv"

# -----------------------------
# Required reference and input files
# -----------------------------
GENOME_FASTA="/path/to/reference/genome.fa"
CHROM_SIZES="/path/to/reference/genome.chrom.sizes"
PEAKS_BED="/path/to/reference/consensus_peaks.bed"
GENES_BED="/path/to/reference/genes.bed"
VARIANT_LIST="/path/to/variants.tsv"

# ChromBPNet variant-scorer source directory and annotation script.
VARIANT_SCORER_SRC="/path/to/variant-scorer/src"
ANNOTATION_SCRIPT="/path/to/variant_annotation.py"

# -----------------------------
# Input schema assumptions
# -----------------------------
# The averaged variant-score TSV is expected to begin with:
#   chromosome, 1-based position, allele1, allele2, variant_id
SCORE_CHROM_COL=1
SCORE_POS_COL=2
SCORE_ALLELE1_COL=3
SCORE_ALLELE2_COL=4
SCORE_ID_COL=5

# FiNeMo hit TSV column numbers used to build the annotation BED.
# Output BED columns are:
#   chromosome, start, end, motif_name, hit_importance, strand, label
HIT_CHROM_COL=1
HIT_START_COL=2
HIT_END_COL=3
HIT_MOTIF_COL=6
HIT_IMPORTANCE_COL=11
HIT_STRAND_COL=13
HIT_LABEL_SUFFIX="modisco_hit"

# -----------------------------
# Runtime controls
# -----------------------------
RUN_SCORING=${RUN_SCORING:-true}
RUN_AVERAGING=${RUN_AVERAGING:-true}
RUN_ANNOTATION=${RUN_ANNOTATION:-true}
RUN_SHAP=${RUN_SHAP:-false}

OVERWRITE_OUTPUTS=${OVERWRITE_OUTPUTS:-false}
REBUILD_ANNOTATION_INPUTS=${REBUILD_ANNOTATION_INPUTS:-true}

SCORING_BATCH_SIZE=${SCORING_BATCH_SIZE:-20000}
SHAP_BATCH_SIZE=${SHAP_BATCH_SIZE:-10000}
SHAP_TYPES=(counts profile)

VARIANT_SCHEMA="chrombpnet"
ANNOTATION_SCHEMA="bed"

################################################################################
# END CONFIGURATION
################################################################################

VARIANT_SCORING_SCRIPT="${VARIANT_SCORER_SRC}/variant_scoring.py"
VARIANT_SUMMARY_SCRIPT="${VARIANT_SCORER_SRC}/variant_summary_across_folds.py"
VARIANT_SHAP_SCRIPT="${VARIANT_SCORER_SRC}/variant_shap.py"

SCORES_DIR="${OUTPUT_DIR}/scores"
AVERAGED_DIR="${OUTPUT_DIR}/averaged_scores"
SHAP_DIR="${OUTPUT_DIR}/variant_shap"

render_template() {
  local template=$1
  local condition=${2:-}
  local fold=${3:-}

  template=${template//\{cell_type\}/${CELL_TYPE}}
  template=${template//\{condition\}/${condition}}
  template=${template//\{fold\}/${fold}}
  printf '%s\n' "${template}"
}

require_file() {
  local file=$1
  local description=$2

  if [[ ! -s "${file}" ]]; then
    echo "ERROR: ${description} is missing or empty: ${file}" >&2
    exit 1
  fi
}

require_dir() {
  local directory=$1
  local description=$2

  if [[ ! -d "${directory}" ]]; then
    echo "ERROR: ${description} does not exist: ${directory}" >&2
    exit 1
  fi
}

is_true() {
  [[ "${1,,}" == "true" ]]
}

initialize_environment() {
  if is_true "${LOAD_MODULES}"; then
    if ! type module >/dev/null 2>&1; then
      echo "ERROR: Environment modules are unavailable, but LOAD_MODULES=true." >&2
      exit 1
    fi

    module purge
    local module_name
    for module_name in "${MODULES[@]}"; do
      module load "${module_name}"
    done
  fi

  if [[ -n "${CONDA_ENV}" ]]; then
    if [[ -n "${CONDA_INIT_SCRIPT}" ]]; then
      require_file "${CONDA_INIT_SCRIPT}" "Conda initialization script"
      # shellcheck source=/dev/null
      source "${CONDA_INIT_SCRIPT}"
    else
      if ! command -v conda >/dev/null 2>&1; then
        echo "ERROR: conda was not found. Load it or set CONDA_INIT_SCRIPT." >&2
        exit 1
      fi
      # shellcheck source=/dev/null
      source "$(conda info --base)/etc/profile.d/conda.sh"
    fi

    conda activate "${CONDA_ENV}"
  fi

  if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi
  fi

  if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    echo "ERROR: Python executable not found: ${PYTHON_BIN}" >&2
    exit 1
  fi
}

validate_configuration() {
  if [[ ${#CONDITIONS[@]} -eq 0 ]]; then
    echo "ERROR: CONDITIONS must contain at least one condition." >&2
    exit 1
  fi

  if [[ ${#FOLDS[@]} -eq 0 ]]; then
    echo "ERROR: FOLDS must contain at least one fold." >&2
    exit 1
  fi

  require_file "${GENOME_FASTA}" "Genome FASTA"
  require_file "${CHROM_SIZES}" "Chromosome sizes file"
  require_file "${PEAKS_BED}" "Peak BED"
  require_file "${VARIANT_LIST}" "Variant list"
  require_dir "${VARIANT_SCORER_SRC}" "Variant-scorer source directory"

  if is_true "${RUN_SCORING}"; then
    require_file "${VARIANT_SCORING_SCRIPT}" "Variant-scoring script"
  fi

  if is_true "${RUN_AVERAGING}"; then
    require_file "${VARIANT_SUMMARY_SCRIPT}" "Fold-summary script"
  fi

  if is_true "${RUN_ANNOTATION}"; then
    require_file "${GENES_BED}" "Gene BED"
    require_file "${ANNOTATION_SCRIPT}" "Variant-annotation script"
    require_dir "${HITS_SOURCE_DIR}" "FiNeMo hit source directory"
  fi

  if is_true "${RUN_SHAP}"; then
    require_file "${VARIANT_SHAP_SCRIPT}" "Variant-SHAP script"
  fi
}

prepare_output_directories() {
  mkdir -p "${SCORES_DIR}" "${AVERAGED_DIR}" "${SHAP_DIR}" \
    "${ANNOTATION_INPUT_DIR}"

  local condition
  for condition in "${CONDITIONS[@]}"; do
    mkdir -p \
      "${AVERAGED_DIR}/${condition}" \
      "${SHAP_DIR}/${condition}"
  done
}

score_variants() {
  echo "========================================"
  echo "STEP 1: Variant scoring per fold"
  echo "========================================"

  local condition fold model_path output_prefix output_file
  for condition in "${CONDITIONS[@]}"; do
    echo "Processing condition: ${condition}"

    for fold in "${FOLDS[@]}"; do
      model_path=$(render_template "${MODEL_PATH_TEMPLATE}" "${condition}" "${fold}")
      output_prefix="${SCORES_DIR}/${condition}_${fold}"
      output_file="${output_prefix}.variant_scores.tsv"

      require_file "${model_path}" "Model for ${condition}, fold ${fold}"

      if [[ -s "${output_file}" ]] && ! is_true "${OVERWRITE_OUTPUTS}"; then
        echo "  Skipping ${condition}, fold ${fold}: output exists"
        continue
      fi

      echo "  Scoring ${condition}, fold ${fold}"
      "${PYTHON_BIN}" "${VARIANT_SCORING_SCRIPT}" \
        --list "${VARIANT_LIST}" \
        --genome "${GENOME_FASTA}" \
        --model "${model_path}" \
        --out_prefix "${output_prefix}" \
        --chrom_sizes "${CHROM_SIZES}" \
        --peaks "${PEAKS_BED}" \
        --peak_genome "${GENOME_FASTA}" \
        -bs "${SCORING_BATCH_SIZE}"
    done
  done
}

average_scores() {
  echo "========================================"
  echo "STEP 2: Average scores across folds"
  echo "========================================"

  local condition fold out_prefix averaged_output score_file
  local -a score_files

  for condition in "${CONDITIONS[@]}"; do
    out_prefix="${AVERAGED_DIR}/${condition}/${condition}_averaged"
    averaged_output="${out_prefix}.mean.variant_scores.tsv"

    if [[ -s "${averaged_output}" ]] && ! is_true "${OVERWRITE_OUTPUTS}"; then
      echo "Skipping ${condition}: averaged output exists"
      continue
    fi

    score_files=()
    for fold in "${FOLDS[@]}"; do
      score_file="${SCORES_DIR}/${condition}_${fold}.variant_scores.tsv"
      require_file "${score_file}" "Variant scores for ${condition}, fold ${fold}"
      score_files+=("${score_file}")
    done

    echo "Averaging scores for ${condition}"
    "${PYTHON_BIN}" "${VARIANT_SUMMARY_SCRIPT}" \
      --score_dir "${SCORES_DIR}" \
      --score_list "${score_files[@]}" \
      --out_prefix "${out_prefix}" \
      --schema "${VARIANT_SCHEMA}"
  done
}

sort_annotation_inputs() {
  local genes_sorted=$1
  local peaks_sorted=$2

  if is_true "${REBUILD_ANNOTATION_INPUTS}" || [[ ! -s "${genes_sorted}" ]]; then
    LC_ALL=C sort -k1,1 -k2,2n "${GENES_BED}" > "${genes_sorted}"
  fi

  if is_true "${REBUILD_ANNOTATION_INPUTS}" || [[ ! -s "${peaks_sorted}" ]]; then
    LC_ALL=C sort -k1,1 -k2,2n "${PEAKS_BED}" > "${peaks_sorted}"
  fi
}

build_hits_union() {
  local condition=$1
  local output_hits=$2
  local hit_glob
  local -a hit_files

  if [[ -s "${output_hits}" ]] && ! is_true "${REBUILD_ANNOTATION_INPUTS}"; then
    echo "Using existing motif-hit union: ${output_hits}"
    return
  fi

  hit_glob=$(render_template "${HIT_FILE_GLOB_TEMPLATE}" "${condition}")
  mapfile -t hit_files < <(compgen -G "${hit_glob}" | LC_ALL=C sort || true)

  if [[ ${#hit_files[@]} -eq 0 ]]; then
    echo "ERROR: No FiNeMo hit files matched for ${condition}:" >&2
    echo "  ${hit_glob}" >&2
    exit 1
  fi

  echo "Building motif-hit union for ${condition} from ${#hit_files[@]} files"

  # Some FiNeMo exports may concatenate a numeric field and a peak-strand token.
  # The sed expression preserves the original pipeline's repair step.
  sed -E 's/([0-9])([+-]peak)/\1\t\2/g' "${hit_files[@]}" | \
    awk \
      -v condition="${condition}" \
      -v label_suffix="${HIT_LABEL_SUFFIX}" \
      -v chrom_col="${HIT_CHROM_COL}" \
      -v start_col="${HIT_START_COL}" \
      -v end_col="${HIT_END_COL}" \
      -v motif_col="${HIT_MOTIF_COL}" \
      -v importance_col="${HIT_IMPORTANCE_COL}" \
      -v strand_col="${HIT_STRAND_COL}" \
      'BEGIN { OFS="\t" }
       $start_col ~ /^[0-9]+$/ && $end_col ~ /^[0-9]+$/ {
         print $chrom_col, $start_col, $end_col, $motif_col, \
               $importance_col, $strand_col, condition "_" label_suffix
       }' | \
    LC_ALL=C sort -k1,1 -k2,2n -u > "${output_hits}"

  require_file "${output_hits}" "Motif-hit union for ${condition}"

  local column_count
  column_count=$(awk -F'\t' 'NR == 1 { print NF; exit }' "${output_hits}")
  if [[ "${column_count}" -ne 7 ]]; then
    echo "ERROR: Expected 7 columns in ${output_hits}; found ${column_count}." >&2
    exit 1
  fi
}

annotate_variants() {
  echo "========================================"
  echo "STEP 3: Annotate averaged variants"
  echo "========================================"

  local genes_sorted="${ANNOTATION_INPUT_DIR}/genes.lexsort.bed"
  local peaks_sorted="${ANNOTATION_INPUT_DIR}/peaks.lexsort.bed"
  sort_annotation_inputs "${genes_sorted}" "${peaks_sorted}"

  local condition hits_bed averaged_scores annotation_list annotation_prefix
  local annotation_output

  for condition in "${CONDITIONS[@]}"; do
    hits_bed="${ANNOTATION_INPUT_DIR}/${condition}.all_folds_hits.bed"
    build_hits_union "${condition}" "${hits_bed}"

    averaged_scores="${AVERAGED_DIR}/${condition}/${condition}_averaged.mean.variant_scores.tsv"
    annotation_list="${AVERAGED_DIR}/${condition}/${condition}_averaged.annotation_input.tsv"
    annotation_prefix="${AVERAGED_DIR}/${condition}/${condition}_averaged.annotated"
    annotation_output="${annotation_prefix}.annotations.tsv"

    require_file "${averaged_scores}" "Averaged scores for ${condition}"

    if [[ -s "${annotation_output}" ]] && ! is_true "${OVERWRITE_OUTPUTS}"; then
      echo "Skipping annotation for ${condition}: output exists"
      continue
    fi

    echo "Preparing annotation input for ${condition}"
    {
      printf 'chr\tpos\tend\tallele1\tallele2\tvariant_id\n'
      awk \
        -v chrom_col="${SCORE_CHROM_COL}" \
        -v pos_col="${SCORE_POS_COL}" \
        -v allele1_col="${SCORE_ALLELE1_COL}" \
        -v allele2_col="${SCORE_ALLELE2_COL}" \
        -v id_col="${SCORE_ID_COL}" \
        'BEGIN { OFS="\t" }
         NR > 1 {
           print $chrom_col, $pos_col - 1, $pos_col, \
                 $allele1_col, $allele2_col, $id_col
         }' "${averaged_scores}" | \
        LC_ALL=C sort -k1,1 -k2,2n
    } > "${annotation_list}"

    require_file "${annotation_list}" "Annotation input for ${condition}"

    if is_true "${OVERWRITE_OUTPUTS}"; then
      rm -f "${annotation_output}"
    fi

    echo "Annotating variants for ${condition}"
    "${PYTHON_BIN}" "${ANNOTATION_SCRIPT}" \
      --list "${annotation_list}" \
      --out_prefix "${annotation_prefix}" \
      --peaks "${peaks_sorted}" \
      --genes "${genes_sorted}" \
      --hits "${hits_bed}" \
      --schema "${ANNOTATION_SCHEMA}"
  done
}

compute_variant_shap() {
  echo "========================================"
  echo "STEP 4: Variant SHAP per fold"
  echo "========================================"

  local condition fold model_path shap_prefix shap_type
  local all_outputs_exist

  for condition in "${CONDITIONS[@]}"; do
    for fold in "${FOLDS[@]}"; do
      model_path=$(render_template "${MODEL_PATH_TEMPLATE}" "${condition}" "${fold}")
      require_file "${model_path}" "Model for ${condition}, fold ${fold}"

      shap_prefix="${SHAP_DIR}/${condition}/${condition}_${fold}"
      all_outputs_exist=true

      for shap_type in "${SHAP_TYPES[@]}"; do
        if [[ ! -s "${shap_prefix}.variant_shap.${shap_type}.h5" ]]; then
          all_outputs_exist=false
          break
        fi
      done

      if is_true "${all_outputs_exist}" && ! is_true "${OVERWRITE_OUTPUTS}"; then
        echo "Skipping SHAP for ${condition}, fold ${fold}: outputs exist"
        continue
      fi

      echo "Computing SHAP for ${condition}, fold ${fold}"
      "${PYTHON_BIN}" "${VARIANT_SHAP_SCRIPT}" \
        --list "${VARIANT_LIST}" \
        --genome "${GENOME_FASTA}" \
        --chrom_sizes "${CHROM_SIZES}" \
        --model "${model_path}" \
        --out_prefix "${shap_prefix}" \
        --schema "${VARIANT_SCHEMA}" \
        --shap_type "${SHAP_TYPES[@]}" \
        --batch_size "${SHAP_BATCH_SIZE}"
    done
  done
}

main() {
  initialize_environment
  validate_configuration
  prepare_output_directories

  echo "Cell type: ${CELL_TYPE}"
  echo "Conditions: ${CONDITIONS[*]}"
  echo "Folds: ${FOLDS[*]}"
  echo "Output directory: ${OUTPUT_DIR}"

  if is_true "${RUN_SCORING}"; then
    score_variants
  else
    echo "Skipping STEP 1: RUN_SCORING=${RUN_SCORING}"
  fi

  if is_true "${RUN_AVERAGING}"; then
    average_scores
  else
    echo "Skipping STEP 2: RUN_AVERAGING=${RUN_AVERAGING}"
  fi

  if is_true "${RUN_ANNOTATION}"; then
    annotate_variants
  else
    echo "Skipping STEP 3: RUN_ANNOTATION=${RUN_ANNOTATION}"
  fi

  if is_true "${RUN_SHAP}"; then
    compute_variant_shap
  else
    echo "Skipping STEP 4: RUN_SHAP=${RUN_SHAP}"
  fi

  echo "Pipeline complete."
}

main "$@"
