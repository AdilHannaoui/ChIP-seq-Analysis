#!/usr/bin/env bash
# ==========================
# RIP-seq IGV Batch Snapshot Module
# Author: Adil Hannaoui Anaaoui
# ==========================

set -euo pipefail

# --------------------------
# Load OUTPUT_DIR from config.R
# --------------------------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_R="${PROJECT_DIR}/config.R"

OUTPUT_DIR=$(Rscript -e "source('${CONFIG_R}'); cat(OUTPUT_DIR)")
SNAP_DIR="${OUTPUT_DIR}/igv_snapshots"
mkdir -p "${SNAP_DIR}"

# --------------------------
# Input files (adjust if needed)
# --------------------------
BAM_WT="${OUTPUT_DIR}/WT_aligned_sorted.bam"
BAM_A="${OUTPUT_DIR}/A_aligned_sorted.bam"
BAM_D="${OUTPUT_DIR}/D_aligned_sorted.bam"

PEAKS_WT="${OUTPUT_DIR}/WT_significant_peaks_with_coords.csv"
PEAKS_A="${OUTPUT_DIR}/A_significant_peaks_with_coords.csv"
PEAKS_D="${OUTPUT_DIR}/D_significant_peaks_with_coords.csv"

GENOME="sacCer3"

# --------------------------
# Generate IGV batch script
# --------------------------
BATCH_FILE="${SNAP_DIR}/igv_batch_script.igv"

echo "new" > "${BATCH_FILE}"
echo "genome ${GENOME}" >> "${BATCH_FILE}"

echo "load ${BAM_WT}" >> "${BATCH_FILE}"
echo "load ${BAM_A}" >> "${BATCH_FILE}"
echo "load ${BAM_D}" >> "${BATCH_FILE}"

echo "snapshotDirectory ${SNAP_DIR}" >> "${BATCH_FILE}"

# --------------------------
# Generate snapshots for each peak
# --------------------------
generate_snapshots () {
    local peaks_file="$1"
    local label="$2"

    tail -n +2 "${peaks_file}" | while IFS=',' read -r chrom start end rest; do
        region="${chrom}:${start}-${end}"
        snapshot_name="${label}_${chrom}_${start}_${end}.png"

        echo "goto ${region}" >> "${BATCH_FILE}"
        echo "snapshot ${snapshot_name}" >> "${BATCH_FILE}"
    done
}

generate_snapshots "${PEAKS_WT}" "WT"
generate_snapshots "${PEAKS_A}"  "A"
generate_snapshots "${PEAKS_D}"  "D"

echo "exit" >> "${BATCH_FILE}"

# --------------------------
# Run IGV in batch mode
# --------------------------
echo "Running IGV batch mode..."
igv.sh -b "${BATCH_FILE}"

echo "IGV snapshots saved in: ${SNAP_DIR}"
