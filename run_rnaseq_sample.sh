#!/usr/bin/env bash
set -uo pipefail

# ============================================================
# RNA-seq per-sample pipeline: trim -> align -> quantify -> multiqc
# Usage: ./run_rnaseq_sample.sh
#   (the script will interactively ask for the R1/R2 fastq
#    files, the GTF annotation file, the STAR genome index
#    directory, and the analysis directory to use)
#
# Non-interactive usage is still supported:
#   ./run_rnaseq_sample.sh <R1.fastq.gz> <R2.fastq.gz> <annotation.gtf[.gz]> <star_index_dir> <analysis_dir>
#
# All outputs (trimmed reads, alignments, counts, multiqc report)
# are centralized under <analysis_dir>/{trimmed,aligned}/<SAMPLE>/
# and <analysis_dir>/multiqc/
# ============================================================

# ---- Paths (edit these once if your layout changes) ----
THREADS_TRIM=8
THREADS_STAR=12
THREADS_FC=8

# ---- Colors for status output (auto-disabled if not a terminal) ----
if [ -t 1 ]; then
    C_GREEN='\033[0;32m'; C_RED='\033[0;31m'; C_YELLOW='\033[1;33m'
    C_BLUE='\033[0;34m'; C_BOLD='\033[1m'; C_RESET='\033[0m'
else
    C_GREEN=''; C_RED=''; C_YELLOW=''; C_BLUE=''; C_BOLD=''; C_RESET=''
fi

CURRENT_STEP="Initialization"
CURRENT_LOG=""
SCRIPT_START_TIME=$(date +%s)

# ---- Error handling ----
# Catches any failing command (since set -e is intentionally NOT used,
# we check exit codes explicitly after each major command instead, which
# lets us print a clean, specific error message rather than a raw trace).
fail() {
    local msg="$1"
    echo ""
    echo -e "${C_RED}${C_BOLD}✗ FAILED at step: ${CURRENT_STEP}${C_RESET}"
    echo -e "${C_RED}  ${msg}${C_RESET}"
    if [ -n "${CURRENT_LOG}" ]; then
        echo -e "${C_RED}  See log file for details: ${CURRENT_LOG}${C_RESET}"
    fi
    exit 1
}

trap 'fail "Unexpected error (exit code $?) near line $LINENO."' ERR

# ---- Helper: print a step header ----
step_header() {
    echo ""
    echo -e "${C_BLUE}${C_BOLD}=== [${SAMPLE:-setup}] $1 ===${C_RESET}"
}

# ---- Helper: format seconds as Hh Mm Ss ----
format_duration() {
    local total_seconds=$1
    printf '%dh %dm %ds' $((total_seconds/3600)) $(((total_seconds%3600)/60)) $((total_seconds%60))
}

# ---- Helper: run a command, streaming its output live while also
#      logging it, showing a spinner-free progress line, and checking
#      its exit code explicitly. ----
run_step() {
    local step_name="$1"; shift
    CURRENT_STEP="${step_name}"
    local step_log="${LOG_DIR}/${step_name// /_}.log"
    CURRENT_LOG="${step_log}"
    local start_ts end_ts duration

    step_header "${step_name}"
    echo "Command: $*" > "${step_log}"
    echo "Started: $(date)" >> "${step_log}"
    echo "---" >> "${step_log}"

    start_ts=$(date +%s)

    # Stream output to both the terminal and the log file in real time.
    set +e
    "$@" 2>&1 | tee -a "${step_log}"
    local exit_code=${PIPESTATUS[0]}
    set -e

    end_ts=$(date +%s)
    duration=$(format_duration $((end_ts - start_ts)))

    if [ "${exit_code}" -ne 0 ]; then
        fail "'${step_name}' exited with code ${exit_code} after ${duration}."
    fi

    echo -e "${C_GREEN}✓ ${step_name} completed in ${duration}${C_RESET}"
}

# ---- Helper: prompt for a file path and validate it exists ----
prompt_for_file() {
    local prompt_text="$1"
    local path=""
    while true; do
        read -e -p "${prompt_text}: " path
        path="${path/#\~/$HOME}"   # expand leading ~
        if [ -f "${path}" ]; then
            echo "${path}"
            return 0
        else
            echo -e "${C_YELLOW}  File not found: '${path}'. Please try again.${C_RESET}" >&2
        fi
    done
}

# ---- Helper: prompt for a directory path and validate it exists ----
prompt_for_dir() {
    local prompt_text="$1"
    local path=""
    while true; do
        read -e -p "${prompt_text}: " path
        path="${path/#\~/$HOME}"   # expand leading ~
        if [ -d "${path}" ]; then
            echo "${path}"
            return 0
        else
            echo -e "${C_YELLOW}  Directory not found: '${path}'. Please try again.${C_RESET}" >&2
        fi
    done
}

# ---- Helper: prompt for an output directory. Creates it if it doesn't
#      exist yet (after confirmation), since analysis dirs are often new. ----
prompt_for_output_dir() {
    local prompt_text="$1"
    local path=""
    local confirm=""
    while true; do
        read -e -p "${prompt_text}: " path
        path="${path/#\~/$HOME}"   # expand leading ~
        if [ -z "${path}" ]; then
            echo -e "${C_YELLOW}  Please enter a path.${C_RESET}" >&2
            continue
        fi
        if [ -d "${path}" ]; then
            echo "${path}"
            return 0
        fi
        read -e -p "  '${path}' does not exist. Create it? [Y/n]: " confirm
        confirm="${confirm:-Y}"
        if [[ "${confirm}" =~ ^[Yy] ]]; then
            if mkdir -p "${path}"; then
                echo "${path}"
                return 0
            else
                echo -e "${C_YELLOW}  Could not create '${path}'. Please try a different path.${C_RESET}" >&2
            fi
        fi
    done
}

# ---- Pre-flight: make sure required tools are installed ----
check_dependencies() {
    local missing=()
    for tool in trim_galore STAR featureCounts multiqc; do
        if ! command -v "${tool}" >/dev/null 2>&1; then
            missing+=("${tool}")
        fi
    done
    if [ ${#missing[@]} -ne 0 ]; then
        echo -e "${C_RED}${C_BOLD}Missing required tool(s): ${missing[*]}${C_RESET}"
        echo "Make sure your conda/mamba environment is activated (e.g. 'conda activate rnaseq_env')."
        exit 1
    fi
}

CURRENT_STEP="Checking dependencies"
check_dependencies

# ---- Get R1, R2, GTF, STAR genome index, and analysis dir either from args or interactively ----
CURRENT_STEP="Collecting inputs"
if [ $# -eq 5 ]; then
    R1="$1"
    R2="$2"
    GTF="$3"
    STAR_INDEX="$4"
    ANALYSIS_DIR="$5"
else
    echo -e "${C_BOLD}=== RNA-seq Pipeline Setup ===${C_RESET}"
    R1=$(prompt_for_file "Path to R1 (forward) fastq.gz file")
    R2=$(prompt_for_file "Path to R2 (reverse) fastq.gz file")
    GTF=$(prompt_for_file "Path to GTF annotation file (.gtf or .gtf.gz)")
    STAR_INDEX=$(prompt_for_dir "Path to STAR genome index directory")
    ANALYSIS_DIR=$(prompt_for_output_dir "Path to analysis directory (outputs will be centralized here)")
fi

[ -f "${R1}" ]  || fail "R1 file not found: ${R1}"
[ -f "${R2}" ]  || fail "R2 file not found: ${R2}"
[ -f "${GTF}" ] || fail "GTF file not found: ${GTF}"
[ -d "${STAR_INDEX}" ] || fail "STAR index directory not found: ${STAR_INDEX}"
[ -n "${ANALYSIS_DIR}" ] || fail "Analysis directory not provided."
mkdir -p "${ANALYSIS_DIR}" || fail "Could not create analysis directory: ${ANALYSIS_DIR}"

# ---- Derive a sample name from the R1 filename ----
# Strips common suffixes like _1.fastq.gz, _R1_001.fastq.gz, etc.
SAMPLE=$(basename "${R1}")
SAMPLE="${SAMPLE%.fastq.gz}"
SAMPLE="${SAMPLE%.fq.gz}"
SAMPLE=$(echo "${SAMPLE}" | sed -E 's/(_R?1)(_001)?$//')

# ---- All outputs live under the centralized analysis directory ----
TRIM_DIR="${ANALYSIS_DIR}/trimmed/${SAMPLE}"
ALIGN_DIR="${ANALYSIS_DIR}/aligned/${SAMPLE}"
MULTIQC_DIR="${ANALYSIS_DIR}/multiqc"
LOG_DIR="${ALIGN_DIR}/logs"

# ---- Make sure output directories exist ----
mkdir -p "${TRIM_DIR}" "${ALIGN_DIR}" "${LOG_DIR}" "${MULTIQC_DIR}"

echo ""
echo -e "${C_BOLD}Sample name detected: ${SAMPLE}${C_RESET}"
echo "R1:            ${R1}"
echo "R2:            ${R2}"
echo "GTF:           ${GTF}"
echo "STAR index:    ${STAR_INDEX}"
echo "Analysis dir:  ${ANALYSIS_DIR}"
echo "Logs:          ${LOG_DIR}"
echo ""

# ============================================================
# Step 1: Trim adapters
# ============================================================
run_step "1_Trimming" trim_galore --fastqc --paired --cores "${THREADS_TRIM}" \
    "${R1}" \
    "${R2}" \
    -o "${TRIM_DIR}"

TRIMMED_R1=$(ls "${TRIM_DIR}"/*_val_1.fq.gz 2>/dev/null | head -n 1)
TRIMMED_R2=$(ls "${TRIM_DIR}"/*_val_2.fq.gz 2>/dev/null | head -n 1)
[ -n "${TRIMMED_R1}" ] || fail "No trimmed R1 file (*_val_1.fq.gz) found in ${TRIM_DIR}"
[ -n "${TRIMMED_R2}" ] || fail "No trimmed R2 file (*_val_2.fq.gz) found in ${TRIM_DIR}"

# ============================================================
# Step 2: Align with STAR
# ============================================================
run_step "2_Alignment" STAR --genomeDir "${STAR_INDEX}" \
    --runThreadN "${THREADS_STAR}" --readFilesIn \
    "${TRIMMED_R1}" \
    "${TRIMMED_R2}" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --readFilesCommand zcat \
    --outFileNamePrefix "${ALIGN_DIR}/"

BAM_FILE="${ALIGN_DIR}/Aligned.sortedByCoord.out.bam"
[ -s "${BAM_FILE}" ] || fail "Expected BAM file not found or empty: ${BAM_FILE}"

# ============================================================
# Step 3: Quantify with featureCounts
# ============================================================
run_step "3_Quantification" featureCounts -T "${THREADS_FC}" -t exon -g gene_name -s 0 \
    -p --countReadPairs \
    -a "${GTF}" \
    -o "${ALIGN_DIR}/featureCounts_exon.txt" \
    "${BAM_FILE}"

[ -s "${ALIGN_DIR}/featureCounts_exon.txt" ] || fail "featureCounts output missing or empty."

# ============================================================
# Step 4: Aggregate QC report with MultiQC
# ============================================================
run_step "4_MultiQC" multiqc "${ANALYSIS_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --force

[ -f "${MULTIQC_DIR}/multiqc_report.html" ] || fail "multiqc_report.html not found in ${MULTIQC_DIR}"

# ============================================================
# Done
# ============================================================
TOTAL_DURATION=$(format_duration $(( $(date +%s) - SCRIPT_START_TIME )))
echo ""
echo -e "${C_GREEN}${C_BOLD}=== [${SAMPLE}] All steps completed successfully in ${TOTAL_DURATION} ===${C_RESET}"
echo "Counts file:    ${ALIGN_DIR}/featureCounts_exon.txt"
echo "MultiQC report: ${MULTIQC_DIR}/multiqc_report.html"
echo "Per-step logs:  ${LOG_DIR}/"