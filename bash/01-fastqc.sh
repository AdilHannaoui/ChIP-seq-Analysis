#!/bin/bash
set -euo pipefail

# Purpose: Perform quality control on raw RIP-seq FASTQ files using FastQC (parallelized)
# Author: Adil Hannaoui Anaaoui

source "$(dirname "$0")/config.sh"

mkdir -p "$OUTPUT_DIR/fastqc_results"
cd "$WORKDIR"

export OUTPUT_DIR FASTQ_DIR

# Function executed by GNU parallel
run_fastqc() {
    FASTQ_FILE="$1"
    SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq)

    echo "Processing sample: $SAMPLE_NAME"

    fastqc "$FASTQ_FILE" \
        -o "$OUTPUT_DIR/fastqc_results" \
        > "$OUTPUT_DIR/${SAMPLE_NAME}.log" 2>&1
}

export -f run_fastqc

# Run in parallel using THREADS from config.sh
parallel -j "$THREADS" run_fastqc ::: "$FASTQ_DIR"/*.fastq
