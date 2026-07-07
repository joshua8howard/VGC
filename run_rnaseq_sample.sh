#!/usr/bin/env bash
# ============================================================
# RNA-seq Sample Pipeline Runner
# ============================================================
# WHAT THIS SCRIPT DOES:
#   1. Checks which bioinformatics tools are installed (from the
#      Dockerfile toolset: fastqc, fastq-screen, multiqc,
#      trim-galore, cutadapt, star, bowtie2, samtools,
#      subread/featureCounts, sra-tools).
#   2. Lets the user pick which pipeline steps to run and which
#      tool to use for each step (only showing tools that are
#      actually installed).
#   3. Lets the user point to their FASTQ files (or SRA IDs).
#   4. Runs the chosen tools, per sample, end to end.
#
# HOW TO RUN:
#   chmod +x run_rnaseq_sample.sh
#   ./run_rnaseq_sample.sh
#
# ENVIRONMENT VARIABLES YOU CAN OVERRIDE:
#   THREADS   Number of CPU threads to use for every tool call.
#             Defaults to 8. Example: THREADS=8 ./run_rnaseq_sample.sh
# ============================================================

# `set -e`         : exit immediately if any command exits non-zero
#                     (stops the pipeline on the first real error).
# `set -u`         : treat use of an unset variable as an error
#                     (catches typos like $THREAD instead of $THREADS).
# `set -o pipefail`: if any command in a pipe (cmd1 | cmd2) fails,
#                     the whole pipe reports failure, not just the
#                     last command. Prevents silently swallowed errors.
set -euo pipefail

# ---------- Colors for terminal output ----------
# These are ANSI escape codes. `echo -e` interprets them to colorize
# text so INFO/WARN/ERROR messages are easy to tell apart at a glance.
GREEN='\033[0;32m'   # used for INFO messages and "found" status
YELLOW='\033[1;33m'  # used for WARN messages
RED='\033[0;31m'     # used for ERROR messages and "not found" status
BOLD='\033[1m'       # used for section headers
NC='\033[0m'         # "No Color" — resets formatting after each message

# ---------- Logging helper functions ----------
# "$*" joins all function arguments into a single space-separated
# string, so log "a" "b" prints "a b".
log()  { echo -e "${GREEN}[INFO]${NC} $*"; }          # normal progress messages
warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }         # non-fatal issues (script keeps going)
err()  { echo -e "${RED}[ERROR]${NC} $*" >&2; }       # fatal issues; sent to stderr (>&2)
hdr()  { echo -e "\n${BOLD}== $* ==${NC}"; }          # section header banner

# ============================================================
# GLOBAL STATE (variables shared across functions)
# ============================================================

# Associative array (like a dictionary): tool name -> 1 (found) or 0 (missing).
# Declared with -A so bash treats it as key/value pairs, not a numeric array.
declare -A TOOL_AVAILABLE

# Associative array: pipeline step key (e.g. "trim") -> chosen tool
# name (e.g. "cutadapt"), or the literal string "skip" if the user
# chose not to run that step.
declare -A STEP_TOOL

# Directory the script is launched from.
WORKDIR="$(pwd)"

# All pipeline output (FASTQC reports, trimmed reads, BAMs, counts,
# MultiQC report) is written under this single folder so nothing
# clutters your working directory.
OUTDIR="${WORKDIR}/rnaseq_output"

# Number of CPU threads passed to every multi-threaded tool call
# (fastqc -t, STAR --runThreadN, bowtie2 -p, samtools -@, etc.).
# "${THREADS:-8}" means: use $THREADS if the user set it in their
# environment, otherwise default to 8.
THREADS="${THREADS:-8}"

# expand_home_path <path>
#   Expands a leading ~ to $HOME while leaving relative and absolute
#   paths unchanged. `read` does not perform shell tilde expansion.
expand_home_path() {
    local path="$1"
    printf '%s' "${path/#\~/$HOME}"
}


# ============================================================
# 1. DEPENDENCY CHECK
#    Figures out which tools are actually installed in the
#    current environment (e.g. inside the Docker container or
#    conda env) before asking the user to choose anything.
# ============================================================

# check_tool <tool_name>
#   Tests whether <tool_name> exists as an executable on the PATH.
#   `command -v` prints the path to the binary if found, or nothing
#   (and a non-zero exit code) if not found. We redirect all output
#   to /dev/null since we only care about the exit code (success/fail).
#   Result is stored in the TOOL_AVAILABLE associative array.
check_tool() {
    local tool="$1"   # first (and only) positional argument to this function
    if command -v "$tool" >/dev/null 2>&1; then
        TOOL_AVAILABLE["$tool"]=1   # 1 = found
    else
        TOOL_AVAILABLE["$tool"]=0   # 0 = not found
    fi
}

# check_all_deps
#   Runs check_tool() for every tool this pipeline knows how to use,
#   then prints a found/not-found summary table for the user to see
#   before making any choices.
check_all_deps() {
    hdr "Checking available tools"

    # List of every binary this script might call. Names must match
    # exactly what's on the PATH (case-sensitive — note STAR is
    # capitalized, matching the actual STAR aligner binary name).
    local tools=(fastqc fastq_screen multiqc trim_galore cutadapt STAR bowtie2 \
                 samtools featureCounts prefetch fasterq-dump)

    # Populate TOOL_AVAILABLE for each tool in the list.
    for t in "${tools[@]}"; do
        check_tool "$t"
    done

    # Print a simple two-column table: tool name | status.
    # printf "%-15s %s\n" left-pads the first column to 15 characters
    # wide so the columns line up neatly regardless of tool name length.
    printf "%-15s %s\n" "TOOL" "STATUS"
    printf "%-15s %s\n" "----" "------"
    for t in "${tools[@]}"; do
        if [[ "${TOOL_AVAILABLE[$t]}" -eq 1 ]]; then
            printf "%-15s ${GREEN}found${NC}\n" "$t"
        else
            printf "%-15s ${RED}not found${NC}\n" "$t"
        fi
    done
    echo
}


# ============================================================
# 2. STEP / TOOL SELECTION MENU
#    Generic helper used for every pipeline stage: shows the user
#    only the tools that are actually installed for that stage,
#    plus a "skip this step" option, and records their choice.
# ============================================================

# select_step_tool <label> <step_key> <candidate_tool_1> [<candidate_tool_2> ...]
#   <label>      : human-readable description shown to the user
#                  (e.g. "Step 4: Adapter / quality trimming")
#   <step_key>   : short internal key used to look up the choice later
#                  (e.g. "trim") — this becomes a key in STEP_TOOL
#   <candidates> : one or more tool names that could fulfil this step
#                  (e.g. trim_galore cutadapt) — the function filters
#                  this list down to only the ones that are installed
select_step_tool() {
    local label="$1"; shift          # grab label, then shift it off $@
    local step_key="$1"; shift       # grab step_key, then shift it off $@
    local candidates=("$@")          # everything remaining = candidate tools

    # Build a filtered list containing only the candidates that were
    # found during check_all_deps(). "${TOOL_AVAILABLE[$tool]:-0}"
    # defaults to 0 if $tool somehow isn't in the array at all.
    local available=()
    for tool in "${candidates[@]}"; do
        if [[ "${TOOL_AVAILABLE[$tool]:-0}" -eq 1 ]]; then
            available+=("$tool")
        fi
    done

    # If none of the candidate tools are installed, there's nothing
    # to choose from — auto-skip this step and warn the user why.
    if [[ ${#available[@]} -eq 0 ]]; then
        warn "No installed tool found for '${label}' — this step will be skipped."
        STEP_TOOL["$step_key"]="skip"
        return
    fi

    # Print a numbered menu: each installed tool, plus a final
    # "skip this step" entry so the user can opt out even if a
    # tool IS available.
    echo -e "\n${BOLD}${label}${NC}"
    local i=1
    local menu=("${available[@]}" "skip this step")
    for opt in "${menu[@]}"; do
        echo "  $i) $opt"
        ((i++))
    done

    # Loop until the user enters a valid number within range.
    # The regex ^[0-9]+$ ensures only digits were typed (rejects
    # blank input, letters, etc. before we try to compare it as a number).
    local choice
    while true; do
        read -rp "Select an option [1-${#menu[@]}]: " choice
        if [[ "$choice" =~ ^[0-9]+$ ]] && (( choice >= 1 && choice <= ${#menu[@]} )); then
            break
        fi
        warn "Invalid choice, try again."
    done

    # Bash arrays are zero-indexed, but the menu shown to the user is
    # 1-indexed, hence the "-1" when looking up the chosen entry.
    local selected="${menu[$((choice-1))]}"
    if [[ "$selected" == "skip this step" ]]; then
        STEP_TOOL["$step_key"]="skip"
    else
        STEP_TOOL["$step_key"]="$selected"
    fi
}


# ============================================================
# 3. FULL PIPELINE CONFIGURATION
#    Calls select_step_tool() once per pipeline stage, in order,
#    then shows a summary and asks for final confirmation.
# ============================================================
configure_pipeline() {
    hdr "Configure Pipeline Steps"
    echo "For each step, pick the tool you want to use (only installed tools are shown)."

    # Each call below defines one pipeline stage. The step_key
    # (2nd argument) is what's used everywhere else in the script
    # to check "did the user enable this step, and with what tool?"

    select_step_tool "Step 1: Download raw reads from SRA (optional)" \
        download prefetch
        # Only candidate: prefetch (the sra-tools download utility).
        # fasterq-dump is used automatically afterward if this step
        # is enabled — it's not a separate user choice.

    select_step_tool "Step 2: Raw read QC" \
        qc fastqc
        # Only candidate: fastqc.

    select_step_tool "Step 3: Contamination / species screening (optional)" \
        screen fastq_screen
        # Only candidate: fastq_screen. Checks reads against reference
        # genomes of common contaminant species.

    select_step_tool "Step 4: Adapter / quality trimming" \
        trim trim_galore cutadapt
        # Two candidates — user picks ONE. trim_galore is a wrapper
        # around cutadapt with automatic adapter detection; cutadapt
        # gives more manual control.

    select_step_tool "Step 5: Alignment to reference genome" \
        align STAR bowtie2
        # Two candidates — user picks ONE. STAR is a splice-aware RNA
        # aligner (standard for RNA-seq); bowtie2 is a general-purpose
        # non-splice-aware aligner.

    select_step_tool "Step 6: Gene-level quantification" \
        quant featureCounts
        # Only candidate: featureCounts (from the subread package).
        # Counts reads per gene using a GTF annotation file.

    select_step_tool "Step 7: Aggregate QC report across all steps" \
        multiqc multiqc
        # Only candidate: multiqc. Scans OUTDIR for logs/reports from
        # every other tool and builds one combined HTML report.

    # samtools sort/index isn't a user-facing choice — it's a
    # mandatory companion step whenever alignment runs (BAM files
    # must be sorted and indexed before featureCounts or viewing in
    # IGV). So we auto-enable it here based on whether alignment
    # was selected and whether samtools is actually installed.
    if [[ "${STEP_TOOL[align]}" != "skip" ]]; then
        if [[ "${TOOL_AVAILABLE[samtools]:-0}" -eq 1 ]]; then
            STEP_TOOL[sort]="samtools"
        else
            warn "samtools not found — BAM sorting/indexing will be skipped even though alignment is enabled."
            STEP_TOOL[sort]="skip"
        fi
    else
        STEP_TOOL[sort]="skip"
    fi

    # Print a final summary table of every step and the tool (or
    # "skip") chosen for it, so the user can review before running.
    hdr "Pipeline Summary"
    printf "%-12s %s\n" "STEP" "TOOL"
    printf "%-12s %s\n" "----" "----"
    for key in download qc screen trim align sort quant multiqc; do
        # "${STEP_TOOL[$key]:-skip}" — defensive default in case a
        # key was somehow never set; shouldn't happen but avoids a
        # crash under `set -u` if it did.
        printf "%-12s %s\n" "$key" "${STEP_TOOL[$key]:-skip}"
    done
    echo

    # Ask for final go-ahead. Pressing Enter with no input defaults
    # to "Y" via the ${confirm:-Y} substitution.
    read -rp "Proceed with this configuration? [Y/n]: " confirm
    confirm="${confirm:-Y}"
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
        err "Aborted by user."
        exit 1
    fi
}


# ============================================================
# 4. FASTQ FILE SELECTION
#    Deliberately runs AFTER tool selection/configuration, so the
#    user commits to a pipeline shape before pointing at data.
# ============================================================

# Per-sample file paths, keyed by sample name.
declare -A SAMPLE_R1              # sample name -> path to R1 (or single-end) FASTQ
declare -A SAMPLE_R2              # sample name -> path to R2 FASTQ, or "" if single-end
SAMPLE_NAMES=()                   # ordered list of sample names (plain array, not associative)

select_fastq_files() {
    hdr "Select FASTQ Input"

    # CASE A: the "download" step was enabled — ask for SRA
    # accession numbers instead of local files. The actual
    # download/convert-to-FASTQ happens later, per sample, inside
    # run_download() during the main pipeline loop.
    if [[ "${STEP_TOOL[download]}" != "skip" ]]; then
        # `-a sra_ids` tells `read` to split the input on whitespace
        # into an array called sra_ids (so multiple accessions can be
        # typed space-separated on one line).
        read -rp "Enter one or more SRA accessions (space-separated): " -a sra_ids
        for acc in "${sra_ids[@]}"; do
            SAMPLE_NAMES+=("$acc")
            # "__SRA__" is a sentinel value (placeholder marker) that
            # run_pipeline() checks for later to know "this sample's
            # FASTQ doesn't exist yet — download it first".
            SAMPLE_R1["$acc"]="__SRA__"
            SAMPLE_R2["$acc"]=""
        done
        log "Will download and process: ${sra_ids[*]}"
        return   # exit the function early; skip the local-directory logic below
    fi

    # CASE B: local FASTQ files. `-e` enables readline (Tab-completion,
    # arrow-key history) while typing the path.
    read -erp "Enter the directory containing your FASTQ files: " fq_dir

    # `read` does NOT perform shell expansions like `~` (home
    # directory shortcut) — that only happens when bash parses a
    # command line directly. So we manually replace a leading ~
    # with the actual $HOME path using parameter expansion:
    #   ${fq_dir/#\~/$HOME}  =  "replace ~ at the START of fq_dir with $HOME"
    fq_dir="$(expand_home_path "$fq_dir")"

    if [[ ! -d "$fq_dir" ]]; then
        err "Directory not found: $fq_dir"
        exit 1
    fi

    # Look for "R1" style paired-end files first.
    # `find ... -maxdepth 1` only looks at the top level of the
    # folder (won't recurse into subfolders).
    # The four -name patterns cover common paired-end naming
    # conventions: *_R1*.fastq.gz / *_R1*.fq.gz / *_1.fastq.gz / *_1.fq.gz
    # `mapfile -t r1_files < <(...)` reads the command's output line
    # by line into the array r1_files (the `-t` strips trailing newlines).
    mapfile -t r1_files < <(find "$fq_dir" -maxdepth 1 -type f \
        \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | sort)

    if [[ ${#r1_files[@]} -eq 0 ]]; then
        # No paired-end pattern matched anything — fall back to
        # treating every FASTQ in the folder as an independent
        # single-end sample.
        mapfile -t all_files < <(find "$fq_dir" -maxdepth 1 -type f \
            \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)

        if [[ ${#all_files[@]} -eq 0 ]]; then
            err "No FASTQ files found in $fq_dir"
            exit 1
        fi

        log "No paired R1/R2 pattern detected — treating all files as single-end."
        for f in "${all_files[@]}"; do
            local base
            base="$(basename "$f")"       # strip directory path, keep filename only
            base="${base%.fastq.gz}"      # strip .fastq.gz extension if present
            base="${base%.fq.gz}"         # strip .fq.gz extension if present
            SAMPLE_NAMES+=("$base")
            SAMPLE_R1["$base"]="$f"
            SAMPLE_R2["$base"]=""         # empty = single-end marker
        done
    else
        # Paired-end files were found — for each R1 file, derive the
        # expected R2 filename and the sample name.
        for r1 in "${r1_files[@]}"; do
            local r2 sample

            # Build the expected R2 path by substituting the R1
            # marker for the R2 marker in the filename. Each of
            # these substitutions targets a different naming style;
            # only the one that actually matches the filename does
            # anything (the others are no-ops).
            r2="${r1/_R1/_R2}"
            r2="${r2/_1.fastq.gz/_2.fastq.gz}"
            r2="${r2/_1.fq.gz/_2.fq.gz}"

            # Derive a clean sample name by stripping everything from
            # the R1 marker onward. "%%" removes the LONGEST matching
            # suffix, so "Sample1_R1_001.fastq.gz" -> "Sample1".
            sample="$(basename "$r1")"
            sample="${sample%%_R1*}"
            sample="${sample%%_1.fastq.gz}"
            sample="${sample%%_1.fq.gz}"

            SAMPLE_NAMES+=("$sample")
            SAMPLE_R1["$sample"]="$r1"

            # Only register R2 if that exact file actually exists on
            # disk — otherwise treat this sample as single-end and
            # warn the user (their R2 might be misnamed or missing).
            if [[ -f "$r2" ]]; then
                SAMPLE_R2["$sample"]="$r2"
            else
                warn "No matching R2 found for $r1 — treating '$sample' as single-end."
                SAMPLE_R2["$sample"]=""
            fi
        done
    fi

    # Print what was detected so the user can sanity-check before
    # the pipeline actually starts running.
    echo
    log "Detected ${#SAMPLE_NAMES[@]} sample(s):"
    for s in "${SAMPLE_NAMES[@]}"; do
        if [[ -n "${SAMPLE_R2[$s]}" ]]; then
            echo "  - $s  (paired-end)"
        else
            echo "  - $s  (single-end)"
        fi
    done
}


# ============================================================
# 5. REFERENCE / ANNOTATION INPUTS
#    Only prompts for the specific reference files each SELECTED
#    tool actually needs — e.g. no GTF prompt if quantification
#    was skipped.
# ============================================================
GENOME_INDEX=""   # path to STAR genome index dir, or bowtie2 index prefix
GTF_FILE=""       # path to GTF gene annotation file (for featureCounts)
SCREEN_CONF=""    # path to fastq_screen config file

collect_reference_inputs() {
    # STAR needs a pre-built genome index DIRECTORY (created ahead of
    # time with `STAR --runMode genomeGenerate`).
    if [[ "${STEP_TOOL[align]}" == "STAR" ]]; then
        read -erp "Path to STAR genome index directory: " GENOME_INDEX
        GENOME_INDEX="$(expand_home_path "$GENOME_INDEX")"

    # bowtie2 needs an index PREFIX (e.g. /path/to/index/genome —
    # bowtie2 appends .1.bt2, .2.bt2, etc. itself), not a directory.
    elif [[ "${STEP_TOOL[align]}" == "bowtie2" ]]; then
        read -erp "Path/prefix to bowtie2 index: " GENOME_INDEX
        GENOME_INDEX="$(expand_home_path "$GENOME_INDEX")"
    fi

    # featureCounts needs a GTF file to know gene/exon coordinates.
    if [[ "${STEP_TOOL[quant]}" == "featureCounts" ]]; then
        read -erp "Path to GTF annotation file: " GTF_FILE
        GTF_FILE="$(expand_home_path "$GTF_FILE")"
    fi

    # fastq_screen needs its own config file listing which reference
    # genomes to screen reads against (paths to bowtie2/bwa indexes
    # for each contaminant species).
    if [[ "${STEP_TOOL[screen]}" == "fastq_screen" ]]; then
        read -erp "Path to fastq_screen config file: " SCREEN_CONF
        SCREEN_CONF="$(expand_home_path "$SCREEN_CONF")"

        if [[ ! -f "$SCREEN_CONF" ]]; then
            err "fastq_screen config not found: $SCREEN_CONF"
            exit 1
        fi
    fi
}


# ============================================================
# 6. PER-SAMPLE STEP RUNNERS
#    Each function below wraps one external tool call. They're
#    invoked once per sample (except run_quant and run_multiqc,
#    which operate across all samples at once).
# ============================================================

# run_download <sra_accession>
#   Downloads raw .sra data via `prefetch`, then converts it to
#   FASTQ with `fasterq-dump`. Updates SAMPLE_R1/SAMPLE_R2 with the
#   resulting FASTQ paths so downstream steps can find them.
run_download() {
    local acc="$1"
    log "[$acc] Downloading from SRA..."

    # prefetch <accession> -O <output_dir>
    #   -O : output directory for the downloaded .sra file
    prefetch "$acc" -O "${OUTDIR}/sra"

    if [[ "${TOOL_AVAILABLE[fasterq-dump]:-0}" -eq 1 ]]; then
        # fasterq-dump <sra_file>
        #   --split-files : write R1/R2 to separate files for
        #                   paired-end data (instead of one interleaved file)
        #   -O            : output directory for resulting FASTQ files
        #   -e            : number of threads to use
        fasterq-dump "${OUTDIR}/sra/${acc}/${acc}.sra" \
            --split-files -O "${OUTDIR}/fastq" -e "$THREADS"

        # fasterq-dump names output files <accession>_1.fastq /
        # <accession>_2.fastq for paired-end data. Register these
        # paths so the rest of the pipeline can use them.
        SAMPLE_R1["$acc"]="${OUTDIR}/fastq/${acc}_1.fastq"
        if [[ -f "${OUTDIR}/fastq/${acc}_2.fastq" ]]; then
            SAMPLE_R2["$acc"]="${OUTDIR}/fastq/${acc}_2.fastq"
        fi
    else
        err "fasterq-dump not available; cannot convert SRA to FASTQ."
        exit 1
    fi
}

# run_qc <sample_name> <r1_path> <r2_path_or_empty>
#   Runs FastQC to generate a quality report for the raw reads.
run_qc() {
    local sample="$1" r1="$2" r2="$3"
    log "[$sample] Running FastQC..."
    mkdir -p "${OUTDIR}/fastqc"

    # fastqc <file(s)>
    #   -t : number of threads (FastQC can process multiple files
    #        in parallel, one thread per file)
    #   -o : output directory for the HTML/zip report
    # If r2 is non-empty (paired-end), pass both files in one call;
    # otherwise just the single-end file.
    if [[ -n "$r2" ]]; then
        fastqc -t "$THREADS" -o "${OUTDIR}/fastqc" "$r1" "$r2"
    else
        fastqc -t "$THREADS" -o "${OUTDIR}/fastqc" "$r1"
    fi
}

# run_screen <sample_name> <r1_path>
#   Runs fastq_screen to check what fraction of reads map to
#   various reference/contaminant genomes (defined in the user's
#   config file). Only run on R1 for speed/simplicity.
run_screen() {
    local sample="$1" r1="$2"
    log "[$sample] Running fastq_screen..."
    mkdir -p "${OUTDIR}/fastq_screen"

    # fastq_screen
    #   --conf    : path to the fastq_screen config file (lists
    #               reference genome index paths to screen against)
    #   --outdir  : where to write the resulting report
    #   --threads : number of threads to use
    fastq_screen --conf "$SCREEN_CONF" --outdir "${OUTDIR}/fastq_screen" --threads "$THREADS" "$r1"
}

# run_trim <sample_name> <r1_path> <r2_path_or_empty>
#   Adapter/quality trims reads using whichever tool the user chose
#   (trim_galore or cutadapt). Sets the global TRIMMED_R1/TRIMMED_R2
#   variables so run_align() knows which files to use next.
run_trim() {
    local sample="$1" r1="$2" r2="$3"
    mkdir -p "${OUTDIR}/trimmed"

    if [[ "${STEP_TOOL[trim]}" == "trim_galore" ]]; then
        log "[$sample] Trimming with Trim Galore..."
        if [[ -n "$r2" ]]; then
            # trim_galore --paired <r1> <r2>
            #   --paired : run in paired-end mode (trims both files
            #              together, keeping read pairs in sync)
            #   --cores  : number of threads
            #   -o       : output directory
            trim_galore --paired --cores "$THREADS" -o "${OUTDIR}/trimmed" "$r1" "$r2"

            # Trim Galore's default output naming convention appends
            # "_val_1"/"_val_2" before the extension for paired-end
            # output files — we reconstruct those expected paths here
            # so downstream steps know where to find the trimmed reads.
            TRIMMED_R1="${OUTDIR}/trimmed/$(basename "${r1%.fastq.gz}")_val_1.fq.gz"
            TRIMMED_R2="${OUTDIR}/trimmed/$(basename "${r2%.fastq.gz}")_val_2.fq.gz"
        else
            trim_galore --cores "$THREADS" -o "${OUTDIR}/trimmed" "$r1"
            # Single-end Trim Galore output is named "_trimmed" instead.
            TRIMMED_R1="${OUTDIR}/trimmed/$(basename "${r1%.fastq.gz}")_trimmed.fq.gz"
            TRIMMED_R2=""
        fi

    elif [[ "${STEP_TOOL[trim]}" == "cutadapt" ]]; then
        log "[$sample] Trimming with cutadapt..."
        TRIMMED_R1="${OUTDIR}/trimmed/${sample}_R1.trimmed.fastq.gz"

        if [[ -n "$r2" ]]; then
            TRIMMED_R2="${OUTDIR}/trimmed/${sample}_R2.trimmed.fastq.gz"
            # cutadapt
            #   -j        : number of threads
            #   -a        : 3' adapter sequence to trim from R1
            #               (AGATCGGAAGAGC = standard Illumina universal adapter)
            #   -A        : 3' adapter sequence to trim from R2
            #               (uppercase -A = "the R2 equivalent of -a")
            #   -o        : output path for trimmed R1
            #   -p        : output path for trimmed R2 (paired mode)
            cutadapt -j "$THREADS" -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
                -o "$TRIMMED_R1" -p "$TRIMMED_R2" "$r1" "$r2"
        else
            TRIMMED_R2=""
            cutadapt -j "$THREADS" -a AGATCGGAAGAGC -o "$TRIMMED_R1" "$r1"
        fi

    else
        # Trimming step was skipped by the user — pass the original,
        # untrimmed files straight through to alignment.
        TRIMMED_R1="$r1"
        TRIMMED_R2="$r2"
    fi
}

# run_align <sample_name> <r1_path> <r2_path_or_empty>
#   Aligns (trimmed) reads to the reference genome using whichever
#   aligner the user chose (STAR or bowtie2). Sets the global
#   ALIGNED_BAM variable with the path to the resulting sorted BAM.
run_align() {
    local sample="$1" r1="$2" r2="$3"
    mkdir -p "${OUTDIR}/aligned/${sample}"

    if [[ "${STEP_TOOL[align]}" == "STAR" ]]; then
        log "[$sample] Aligning with STAR..."

        # STAR
        #   --runThreadN       : number of threads
        #   --genomeDir        : path to the pre-built STAR genome index
        #   --readFilesIn      : input FASTQ file(s) — one for
        #                        single-end, two (R1 R2) for paired-end
        #   --readFilesCommand : command STAR uses to decompress input
        #                        on the fly (zcat, since inputs are .gz)
        #   --outSAMtype       : output format — BAM, already
        #                        coordinate-sorted (skips a separate
        #                        samtools sort step for STAR output)
        #   --outFileNamePrefix: prefix for all output file names
        #                        (STAR appends fixed suffixes like
        #                        "Aligned.sortedByCoord.out.bam")
        if [[ -n "$r2" ]]; then
            STAR --runThreadN "$THREADS" --genomeDir "$GENOME_INDEX" \
                --readFilesIn "$r1" "$r2" --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix "${OUTDIR}/aligned/${sample}/"
        else
            STAR --runThreadN "$THREADS" --genomeDir "$GENOME_INDEX" \
                --readFilesIn "$r1" --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix "${OUTDIR}/aligned/${sample}/"
        fi
        # STAR's fixed output filename when --outSAMtype BAM SortedByCoordinate is used.
        ALIGNED_BAM="${OUTDIR}/aligned/${sample}/Aligned.sortedByCoord.out.bam"

    elif [[ "${STEP_TOOL[align]}" == "bowtie2" ]]; then
        log "[$sample] Aligning with bowtie2..."
        local sam="${OUTDIR}/aligned/${sample}/${sample}.sam"

        # bowtie2
        #   -p : number of threads
        #   -x : path/prefix to the bowtie2 index
        #   -1 / -2 : R1/R2 input files (paired-end mode)
        #   -U : single input file (single-end mode)
        #   -S : output SAM file path
        # Unlike STAR, bowtie2 outputs plain SAM (not sorted BAM), so
        # we manually sort it into BAM format with samtools afterward.
        if [[ -n "$r2" ]]; then
            bowtie2 -p "$THREADS" -x "$GENOME_INDEX" -1 "$r1" -2 "$r2" -S "$sam"
        else
            bowtie2 -p "$THREADS" -x "$GENOME_INDEX" -U "$r1" -S "$sam"
        fi

        ALIGNED_BAM="${OUTDIR}/aligned/${sample}/${sample}.sorted.bam"

        # samtools sort
        #   -@ : number of threads
        #   -o : output BAM path
        #   <sam> : input file to sort (converts SAM -> sorted BAM in one step)
        samtools sort -@ "$THREADS" -o "$ALIGNED_BAM" "$sam"
        rm -f "$sam"   # delete the intermediate uncompressed SAM to save disk space
    fi
}

# run_sort_index <sample_name>
#   Indexes the sorted BAM produced by run_align(). A BAM index
#   (.bai file) is required by featureCounts, IGV, and most other
#   downstream BAM-consuming tools for fast random access.
run_sort_index() {
    local sample="$1"
    if [[ "${STEP_TOOL[sort]}" == "samtools" && -n "${ALIGNED_BAM:-}" ]]; then
        log "[$sample] Indexing BAM with samtools..."
        # samtools index <bam> — creates a <bam>.bai file alongside it.
        samtools index "$ALIGNED_BAM"
    fi
}

# run_quant
#   Runs featureCounts once per sample so each aligned BAM gets its
#   own exon-level count file alongside the alignment output.
run_quant() {
    log "Running featureCounts across all samples..."
    mkdir -p "${OUTDIR}/counts"

    for s in "${SAMPLE_NAMES[@]}"; do
        local bam_files=()
        local bam

        while IFS= read -r bam; do
            [[ -n "$bam" ]] && bam_files+=("$bam")
        done < <(find "${OUTDIR}/aligned/${s}" -maxdepth 1 -type f -name "*.bam" | sort)

        if [[ ${#bam_files[@]} -eq 0 ]]; then
            warn "[$s] No BAM files found for quantification — skipping featureCounts."
            continue
        fi

        if [[ ${#bam_files[@]} -gt 1 ]]; then
            warn "[$s] Multiple BAM files found; using ${bam_files[0]}"
        fi

        bam="${bam_files[0]}"

        local out_file="${OUTDIR}/counts/${s}_featureCounts_exon.txt"

        if [[ -n "${SAMPLE_R2[$s]:-}" ]]; then
            log "[$s] Counting paired-end reads with featureCounts..."
            # featureCounts
            #   -T : number of threads
            #   -t : feature type to count (exon)
            #   -g : attribute to group by in the GTF (gene_name)
            #   -s : strand-specificity mode (0 = unstranded)
            #   -p : count paired-end fragments instead of single reads
            #   -B : require both ends of each pair to be properly aligned
            #   -C : ignore chimeric fragments with mates on different chromosomes
            #   -a : path to the GTF annotation file
            #   -o : output file path for the count table
            featureCounts -T "$THREADS" -t exon -g gene_name -s 0 -p -B -C \
                -a "$GTF_FILE" \
                -o "$out_file" \
                "$bam"
        else
            log "[$s] Counting single-end reads with featureCounts..."
            # featureCounts
            #   -T : number of threads
            #   -t : feature type to count (exon)
            #   -g : attribute to group by in the GTF (gene_name)
            #   -s : strand-specificity mode (0 = unstranded)
            #   -a : path to the GTF annotation file
            #   -o : output file path for the count table
            featureCounts -T "$THREADS" -t exon -g gene_name -s 0 \
                -a "$GTF_FILE" \
                -o "$out_file" \
                "$bam"
        fi
    done
}

# run_multiqc
#   Scans the entire output directory for log/report files produced
#   by every other tool (FastQC, fastq_screen, Trim Galore/cutadapt,
#   STAR/bowtie2, featureCounts) and aggregates them into a single
#   interactive HTML report.
run_multiqc() {
    log "Aggregating QC reports with MultiQC..."
    # multiqc <search_dir> -o <output_dir>
    #   Recursively scans <search_dir> for recognizable log files
    #   and writes a combined "multiqc_report.html" to -o.
    multiqc "$OUTDIR" -o "$OUTDIR"
}


# ============================================================
# 7. MAIN PIPELINE EXECUTION
#    Loops over every sample and runs the enabled steps in order,
#    passing each step's output as the next step's input. Steps
#    the user chose to skip are simply not called.
# ============================================================
run_pipeline() {
    mkdir -p "$OUTDIR"

    for sample in "${SAMPLE_NAMES[@]}"; do
        hdr "Processing sample: $sample"

        # If this sample is an SRA accession awaiting download
        # (marked with the "__SRA__" sentinel earlier), download and
        # convert it to FASTQ now, before anything else can run.
        if [[ "${STEP_TOOL[download]}" != "skip" && "${SAMPLE_R1[$sample]}" == "__SRA__" ]]; then
            run_download "$sample"
        fi

        # Current best-known R1/R2 paths for this sample (updated by
        # run_download if applicable).
        local r1="${SAMPLE_R1[$sample]}"
        local r2="${SAMPLE_R2[$sample]:-}"

        # Each line below only calls its step function if the user
        # didn't set that step to "skip" during configure_pipeline().
        [[ "${STEP_TOOL[qc]}" != "skip" ]] && run_qc "$sample" "$r1" "$r2"
        [[ "${STEP_TOOL[screen]}" != "skip" ]] && run_screen "$sample" "$r1"

        # Default: if trimming is skipped, the "trimmed" files are
        # just the original files, so alignment still has valid input.
        TRIMMED_R1="$r1"
        TRIMMED_R2="$r2"
        [[ "${STEP_TOOL[trim]}" != "skip" ]] && run_trim "$sample" "$r1" "$r2"

        # Alignment (and its mandatory sort/index companion) only
        # runs if the user enabled it. ALIGNED_BAM is reset to empty
        # each sample so a stale path from a previous sample can't
        # leak into run_sort_index() by mistake.
        ALIGNED_BAM=""
        if [[ "${STEP_TOOL[align]}" != "skip" ]]; then
            run_align "$sample" "$TRIMMED_R1" "$TRIMMED_R2"
            run_sort_index "$sample"
        fi
    done

    # Quantification and MultiQC operate across ALL samples at once,
    # so they run once after the per-sample loop finishes — not
    # inside it.
    [[ "${STEP_TOOL[quant]}" != "skip" ]] && run_quant
    [[ "${STEP_TOOL[multiqc]}" != "skip" ]] && run_multiqc

    hdr "Pipeline complete"
    log "Results saved to: $OUTDIR"
}


# ============================================================
# ENTRY POINT
#    This is the actual order of operations when the script runs.
# ============================================================
main() {
    hdr "RNA-seq Sample Pipeline"
    check_all_deps           # Step 1: see what's installed
    configure_pipeline       # Step 2: user picks tools per step
    select_fastq_files       # Step 3: user points at FASTQ data (or SRA IDs)
    collect_reference_inputs # Step 4: gather any reference files the chosen tools need
    run_pipeline             # Step 5: actually run everything
}

# "$@" forwards any command-line arguments given to the script on
# to main() — currently unused, but keeps the door open for future
# flags (e.g. --threads 8) without changing this line.
main "$@"