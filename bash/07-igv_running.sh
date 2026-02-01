#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq IGV Batch Snapshot Module (trimmed-only, max 3 reps)
# Fixed: IP and IN together, coverage-only, with diagnostics
# ==========================

source "$(dirname "$0")/config.sh"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_R="R/config.R"

OUTPUT_DIR=$(Rscript -e "source('${CONFIG_R}'); cat(OUTPUT_DIR)")
BAM_DIR="${OUTPUT_DIR}/bowtie2"
if [[ ! -d "${BAM_DIR}" ]]; then
    BAM_DIR="${OUTPUT_DIR}"
fi

SNAP_DIR="${OUTPUT_DIR}/igv_snapshots"
mkdir -p "${SNAP_DIR}"
mkdir -p "${OUTPUT_DIR}/logs"

# THREADS and TOP_N (override via env if desired)
THREADS="${THREADS:-4}"
TOP_N="${TOP_N:-10}"

# --------------------------
# Inputs / patterns
# --------------------------
GENOME="sacCer3"
PEAKS_WT="${OUTPUT_DIR}/WT_significant_peaks_with_coords.csv"
PEAKS_M1="${OUTPUT_DIR}/M1_significant_peaks_with_coords.csv"
PEAKS_M12="${OUTPUT_DIR}/M12_significant_peaks_with_coords.csv"

# --------------------------
# Utilities and early files
# --------------------------
log() { echo "[$(date +'%F %T')] $*"; }
err() { echo "ERROR: $*" >&2; }

BATCH_FILE="${SNAP_DIR}/igv_batch_script.igv"
: > "${BATCH_FILE}" || { err "Cannot create/truncate ${BATCH_FILE}"; exit 1; }

# Detect IGV and xvfb
IGV_CMD=""
if command -v igv.sh >/dev/null 2>&1; then
    IGV_CMD="$(command -v igv.sh)"
elif command -v igv >/dev/null 2>&1; then
    IGV_CMD="$(command -v igv)"
elif [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/igv.sh" ]]; then
    IGV_CMD="${CONDA_PREFIX}/bin/igv.sh"
fi

if [[ -z "${IGV_CMD}" ]]; then
    err "igv.sh or igv not found in PATH or CONDA_PREFIX. Adjust PATH or install IGV."
    exit 1
fi

XVFB_CMD=""
if [[ -z "${DISPLAY:-}" ]]; then
    if command -v xvfb-run >/dev/null 2>&1; then
        XVFB_CMD="xvfb-run -a"
        log "No DISPLAY detected: using xvfb-run to execute IGV in headless mode."
    else
        log "No DISPLAY detected and xvfb-run not available. IGV may fail if X is required."
    fi
fi

for cmd in samtools parallel; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        err "Required command '${cmd}' not found in PATH."
        exit 1
    fi
done

# --------------------------
# Helper: normalize chromosome name
# --------------------------
normalize_chr() {
    local chr="$1"
    chr=$(echo "$chr" | sed 's/^[[:space:]"]*//;s/[[:space:]"]*$//')
    if [[ "$chr" =~ ^chr ]]; then
        echo "$chr"
        return
    fi
    echo "chr${chr}"
}

# --------------------------
# Helper: calculate max coverage in a region for a BAM file
# --------------------------
get_max_coverage_in_region() {
    local bam="$1"
    local region="$2"
    
    if [[ ! -f "$bam" ]]; then
        echo "0"
        return
    fi
    
    local max_cov=$(samtools depth -r "$region" "$bam" 2>/dev/null | awk 'BEGIN{max=0} {if($3>max)max=$3} END{print max}')
    
    if [[ -z "$max_cov" || "$max_cov" == "0" ]]; then
        echo "0"
    else
        echo "$max_cov"
    fi
}

# --------------------------
# Merge helper
# --------------------------
merge_and_index() {
    local out="$1"; shift
    local inputs=("$@")
    local exist=()
    for f in "${inputs[@]}"; do [[ -f "$f" ]] && exist+=("$f"); done

    if [[ ${#exist[@]} -eq 0 ]]; then
        log "No BAMs found for ${out}, skipping."
        return 1
    fi

    log "Merging BAMs into ${out}:"
    for f in "${exist[@]}"; do log "  - ${f}"; done

    if [[ ${#exist[@]} -eq 1 ]]; then
        local src="${exist[0]}"
        log "Only one BAM for ${out}: ${src} -> sorting and indexing"
        samtools view -H "${src}" >/dev/null 2>&1 || { err "Invalid BAM: ${src}"; return 1; }
        if [[ "$(realpath "${src}")" == "$(realpath "${out}")" ]]; then
            samtools index "${out}"
            return 0
        fi
        samtools sort -@ "${THREADS}" -o "${out}" "${src}"
        samtools index "${out}"
        return 0
    fi

    log "Merging ${#exist[@]} BAMs -> ${out}"
    local tmp_unsorted="${out}.unsorted.bam"
    samtools merge -@ "${THREADS}" -o "${tmp_unsorted}" "${exist[@]}"
    samtools sort -@ "${THREADS}" -o "${out}" "${tmp_unsorted}"
    samtools index "${out}"
    rm -f "${tmp_unsorted}"
    log "Merge completed: ${out}"
    return 0
}

# --------------------------
# find_bams_for_trimmed: only *_trimmed.bam, fallback to any *.bam
# returns up to MAX_REPS per type (default 3)
# --------------------------
MAX_REPS=3

find_bams_for_trimmed() {
    local cond="$1"
    local type="$2"
    local lc_type
    lc_type="$(echo "$type" | tr '[:upper:]' '[:lower:]')"

    shopt -s nullglob
    local files=()
    local all_files=()

    for f in "${BAM_DIR}/${cond}"*"_trimmed.bam"; do
        [[ -f "$f" ]] && all_files+=("$f")
    done

    if [[ ${#all_files[@]} -eq 0 ]]; then
        for p in "${BAM_DIR}/${cond}"*"_${type}"*".bam" "${BAM_DIR}/${cond}"*"${type}"*".bam" "${BAM_DIR}/${cond}"*"_${lc_type}"*".bam"; do
            for f in $p; do [[ -f "$f" ]] && all_files+=("$f"); done
        done
    fi

    if [[ ${#all_files[@]} -eq 0 ]]; then
        for f in "${BAM_DIR}/${cond}"*".bam"; do [[ -f "$f" ]] && all_files+=("$f"); done
    fi

    for f in "${all_files[@]}"; do
        base="$(basename "$f")"
        if [[ "$type" == "IP" ]]; then
            if [[ "$base" =~ _[Ii][Pp]([0-9]|_|\.|$) ]]; then files+=("$f"); fi
        else
            if [[ "$base" =~ _[Ii][Nn]([0-9]|_|\.|$) ]] || [[ "$base" =~ [Ii][Nn][Pp][Uu][Tt] ]]; then files+=("$f"); fi
        fi
    done

    if [[ ${#files[@]} -gt 0 ]]; then
        mapfile -t files < <(printf "%s\n" "${files[@]}" | awk '!seen[$0]++' | sort)
    fi

    if [[ ${#files[@]} -gt ${MAX_REPS} ]]; then
        files=("${files[@]:0:${MAX_REPS}}")
    fi

    shopt -u nullglob
    printf "%s\n" "${files[@]}"
}

# --------------------------
# Top peaks helper - sort by pvalue and take top N
# --------------------------
get_top_peaks() {
    local infile="$1"
    local n="${2:-$TOP_N}"
    local outfile="$3"
    
    if [[ ! -f "$infile" ]]; then
        echo "WARNING: get_top_peaks: input not found: $infile" >&2
        return 1
    fi
    
    header="$(head -n1 "$infile" || true)"
    if [[ -z "$header" ]]; then
        echo "WARNING: empty header in $infile" >&2
        return 1
    fi
    
    pvalue_idx=$(echo "$header" | awk -F',' '{
        for(i=1; i<=NF; i++) {
            h=$i;
            gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", h);
            lh=tolower(h);
            if(lh=="pvalue" || lh=="p_value" || lh=="p.value" || lh=="pval") {
                print i;
                exit;
            }
        }
    }')
    
    if [[ -z "$pvalue_idx" ]]; then
        log "WARNING: No pvalue column found in $infile, using first $n rows"
        { echo "$header"; tail -n +2 "$infile" | awk 'NF' | head -n "$n"; } > "$outfile"
    else
        log "Sorting by pvalue (column $pvalue_idx) and taking top $n peaks"
        { 
            echo "$header"
            tail -n +2 "$infile" | awk 'NF' | sort -t',' -k"${pvalue_idx}","${pvalue_idx}" -g | head -n "$n"
        } > "$outfile"
    fi
    
    return 0
}

# --------------------------
# Prepare merges per condition
# --------------------------
conds=(WT M1 M12)
declare -A MERGED_IP MERGED_IN

for c in "${conds[@]}"; do
    mapfile -t ips < <(find_bams_for_trimmed "${c}" "IP" || true)
    mapfile -t ins < <(find_bams_for_trimmed "${c}" "IN" || true)

    log "Detected for ${c} IP: ${ips[*]:-<none>}"
    log "Detected for ${c} IN: ${ins[*]:-<none>}"

    out_ip="${OUTPUT_DIR}/${c}_IP_merged_sorted.bam"
    out_in="${OUTPUT_DIR}/${c}_IN_merged_sorted.bam"

    if [[ ${#ips[@]} -gt 0 ]]; then
        merge_and_index "${out_ip}" "${ips[@]}" || log "IP merge failed for ${c}"
        MERGED_IP["${c}"]="${out_ip}"
    fi

    if [[ ${#ins[@]} -gt 0 ]]; then
        merge_and_index "${out_in}" "${ins[@]}" || log "IN merge failed for ${c}"
        MERGED_IN["${c}"]="${out_in}"
    fi
done

# --------------------------
# CRITICAL DIAGNOSTIC: Verify IP and IN are different
# --------------------------
log "=== CRITICAL DIAGNOSTIC: Checking IP vs IN differences ==="
for c in "${conds[@]}"; do
    ip_bam="${MERGED_IP[$c]:-}"
    in_bam="${MERGED_IN[$c]:-}"
    
    if [[ -n "${ip_bam}" && -f "${ip_bam}" && -n "${in_bam}" && -f "${in_bam}" ]]; then
        
        if [[ "$(realpath "${ip_bam}")" == "$(realpath "${in_bam}")" ]]; then
            err "ERROR! IP and IN point to the SAME file for ${c}:"
            err "  IP: ${ip_bam}"
            err "  IN: ${in_bam}"
            continue
        fi
        
        ip_size=$(stat -f%z "${ip_bam}" 2>/dev/null || stat -c%s "${ip_bam}" 2>/dev/null || echo "0")
        in_size=$(stat -f%z "${in_bam}" 2>/dev/null || stat -c%s "${in_bam}" 2>/dev/null || echo "0")
        
        log "${c} - File sizes:"
        log "  IP: ${ip_bam} (${ip_size} bytes)"
        log "  IN: ${in_bam} (${in_size} bytes)"
        
        ip_count=$(samtools view -c "${ip_bam}" 2>/dev/null || echo "0")
        in_count=$(samtools view -c "${in_bam}" 2>/dev/null || echo "0")
        
        log "${c} - Total reads:"
        log "  IP: ${ip_count} reads"
        log "  IN: ${in_count} reads"
        
        if [[ "$ip_count" == "$in_count" && "$ip_count" != "0" ]]; then
            log "⚠️ WARNING: IP and IN have the SAME number of reads for ${c}!"
        fi
        
        log "Comparing first reads for ${c}..."
        ip_first=$(samtools view "${ip_bam}" 2>/dev/null | head -1 | cut -f1 || echo "NONE")
        in_first=$(samtools view "${in_bam}" 2>/dev/null | head -1 | cut -f1 || echo "NONE")
        
        if [[ "$ip_first" == "$in_first" && "$ip_first" != "NONE" ]]; then
            log "⚠️ First reads have the same name for ${c}"
        else
            log "✓ IP and IN appear different for ${c}"
        fi
    fi
done

# --------------------------
# Build batch: IP and IN together, coverage only
# --------------------------
declare -A PEAKS_MAP
PEAKS_MAP[WT]="${PEAKS_WT}"
PEAKS_MAP[M1]="${PEAKS_M1}"
PEAKS_MAP[M12]="${PEAKS_M12}"

tmp_dir="${SNAP_DIR}/tmp_top_peaks"
mkdir -p "${tmp_dir}"

{
    echo "new"
    echo "genome ${GENOME}"
    echo "snapshotDirectory ${SNAP_DIR}"
    echo "maxPanelHeight 1000"
    echo "preference SAM.SHOW_COV_TRACK true"
    echo "preference SAM.SHOW_ALIGNMENT_TRACK false"
    echo "preference SAM.AUTOSCALE true"
} > "${BATCH_FILE}"

for c in "${conds[@]}"; do
    ip_bam="${MERGED_IP[$c]:-}"
    in_bam="${MERGED_IN[$c]:-}"

    peaks_csv="${PEAKS_MAP[$c]}"
    if [[ ! -f "${peaks_csv}" ]]; then
        log "Peak CSV not found for ${c}: ${peaks_csv}, skipping."
        continue
    fi

    topf="${tmp_dir}/${c}_top${TOP_N}.csv"
    get_top_peaks "${peaks_csv}" "${TOP_N}" "${topf}" || { 
        log "Failed to extract top peaks for ${c}, skipping"
        continue
    }
    
    log "Top ${TOP_N} peaks for ${c} saved in: ${topf}"

    tracks_loaded=0
    if [[ -n "${ip_bam}" && -f "${ip_bam}" ]]; then
        echo "load ${ip_bam}" >> "${BATCH_FILE}"
        tracks_loaded=$((tracks_loaded + 1))
        log "Loading IP: ${ip_bam}"
    else
        log "⚠️ Warning: IP merge not found for ${c}"
    fi
    
    if [[ -n "${in_bam}" && -f "${in_bam}" ]]; then
        echo "load ${in_bam}" >> "${BATCH_FILE}"
        tracks_loaded=$((tracks_loaded + 1))
        log "Loading IN: ${in_bam}"
    else
        log "⚠️ Warning: IN merge not found for ${c}"
    fi
    
    if [[ $tracks_loaded -eq 0 ]]; then
        log "No tracks to load for ${c}, skipping"
        continue
    fi
    
    echo "viewaspairs" >> "${BATCH_FILE}"
    echo "collapse" >> "${BATCH_FILE}"

    tail -n +2 "${topf}" | while IFS= read -r line; do
        [[ -z "$line" ]] && continue
        
        chrom=$(echo "$line" | awk -F',' '{val=$2; gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", val); print val;}')
        start=$(echo "$line" | awk -F',' '{val=$3; gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", val); print val;}')
        end=$(echo "$line" | awk -F',' '{val=$4; gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", val); print val;}')
        
        [[ -z "${chrom}" || -z "${start}" || -z "${end}" ]] && continue
        
        if ! [[ "$start" =~ ^[0-9]+$ ]] || ! [[ "$end" =~ ^[0-9]+$ ]]; then
            continue
        fi
        
        chrom_norm=$(normalize_chr "$chrom")
        region="${chrom_norm}:${start}-${end}"
        snapshot="${c}_${chrom}_${start}_${end}.png"
        
        max_ip=0
        max_in=0
        
        if [[ -n "${ip_bam}" && -f "${ip_bam}" ]]; then
            max_ip=$(get_max_coverage_in_region "${ip_bam}" "${region}")
        fi
        
        if [[ -n "${in_bam}" && -f "${in_bam}" ]]; then
            max_in=$(get_max_coverage_in_region "${in_bam}" "${region}")
        fi
        
        max_cov=$((max_ip > max_in ? max_ip : max_in))
        
        max_scale=$(awk "BEGIN {print int($max_cov * 1.1) + 1}")
        
        log "Snapshot: ${region} -> ${snapshot} (IP_max=${
