#!/bin/bash
set -euo pipefail

# Purpose: Perform quality control on raw RIP-seq FASTQ files using FastQC
# Author: Adil Hannaoui Anaaoui

WORKDIR="/mnt/c/Users/rip-seq/"
FASTQ_DIR="$WORKDIR/data/"
OUTPUT_DIR="$WORKDIR/output"
THREADS=6

mkdir -p "$OUTPUT_DIR/fastqc_results"
cd "$WORKDIR"

for FASTQ_FILE in "$FASTQ_DIR"/*.fastq; do
    SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq)
    echo "Processing sample: $SAMPLE_NAME"

    if ! fastqc "$FASTQ_FILE" -o "$OUTPUT_DIR/fastqc_results" > "$OUTPUT_DIR/$SAMPLE_NAME.log" 2>&1; then
        echo "Error: FastQC failed for $SAMPLE_NAME"
        continue
    fi
done
