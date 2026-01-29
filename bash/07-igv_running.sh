#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq IGV Batch Snapshot Module (parallel batch generation)
# Author: Adil Hannaoui Anaaoui
# ==========================

source "$(dirname "$0")/config.sh"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_R="R/config.R"

OUTPUT_DIR=$(Rscript -e "source('${CONFIG_R}'); cat(OUTPUT_DIR)")
SNAP_DIR="${OUTPUT_DIR}/igv_snapshots"
mkdir -p "${SNAP_DIR}"
mkdir -p "${OUTPUT_DIR}/logs"

# --------------------------
# Input files
# --------------------------
BAM_WT="${OUTPUT_DIR}/WT_aligned_sorted.bam"
BAM_A="${OUTPUT_DIR}/A_aligned_sorted.bam"
BAM_D="${OUTPUT_DIR}/D_aligned_sorted.bam"

PEAKS_WT="${OUTPUT_DIR}/WT_significant_peaks_with_coords.csv"
PEAKS_A="${OUTPUT_DIR}/A_significant_peaks_with_coords.csv"
PEAKS_D="${OUTPUT_DIR}/D_significant_peaks_with_coords.csv"

GENOME="sacCer3"

# --------------------------
# Batch script
# --------------------------
BATCH_FILE="${SNAP_DIR}/igv_batch_script.igv"

{
    echo "new"
    echo "genome ${GENOME}"
    echo "load ${BAM_WT}"
    echo "load ${BAM_A}"
    echo "load ${BAM_D}"
    echo "snapshotDirectory ${SNAP_DIR}"
} > "${BATCH_FILE}"

# --------------------------
# Function to generate snapshot commands
# --------------------------
generate_snapshot_lines() {
    local peaks_file="$1"
    local label="$2"

    tail -n +2 "$peaks_file" | while IFS=',' read -r chrom start end rest; do
        region="${chrom}:${start}-${end}"
        snapshot="${label}_${chrom}_${start}_${end}.png"
        echo "goto ${region}"
        echo "snapshot ${snapshot}"
    done
}

export -f generate_snapshot_lines

# --------------------------
# Generate snapshot commands in parallel
# --------------------------
parallel -j "$THREADS" generate_snapshot_lines ::: \
    "$PEAKS_WT" "$PEAKS_A" "$PEAKS_D" ::: \
    "WT" "A" "D" >> "${BATCH_FILE}"

echo "exit" >> "${BATCH_FILE}"

# --------------------------
# Run IGV batch mode
# --------------------------
echo "Running IGV batch mode..."
igv.sh -b "${BATCH_FILE}" > "${OUTPUT_DIR}/logs/igv_batch.log" 2>&1

echo "IGV snapshots saved in: ${SNAP_DIR}"
