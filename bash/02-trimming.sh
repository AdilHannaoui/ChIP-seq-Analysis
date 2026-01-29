#!/bin/bash

# Purpose: Perform adapter trimming and quality filtering of raw RIP-seq FASTQ files
# using Trimmomatic. The output are cleaned FASTQ files suitable for downstream
# alignment.

# Author: Adil Hannaoui Anaaoui

set -euo pipefail

source "$(dirname "$0")/config.sh"

cd "$WORKDIR"
mkdir -p "$OUTPUT_DIR/fastqc_trimmed"

for FASTQ_FILE in "$FASTQ_DIR"/*.fastq; do
    SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq)
    echo "Processing Sample: $SAMPLE_NAME"

    TRIMMED_FASTQ_FILE="$FASTQ_TRIM/${SAMPLE_NAME}_trimmed.fastq"

    # Run Trimmomatic
    java -jar "$TRIMMO_JAR" SE -threads "$THREADS" \
        "$FASTQ_FILE" "$TRIMMED_FASTQ_FILE" \
        ILLUMINACLIP:"$WORKDIR/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10" \
        SLIDINGWINDOW:4:20 MINLEN:20 -phred33

    # Check if trimming was successful
    if [ ! -s "$TRIMMED_FASTQ_FILE" ]; then
        echo "Error: Trimmomatic failed for $SAMPLE_NAME"
        exit 1
    fi

    # FastQC post-trimming
    fastqc "$TRIMMED_FASTQ_FILE" -o "$OUTPUT_DIR/fastqc_trimmed" \
        > "$OUTPUT_DIR/${SAMPLE_NAME}_trimmed.log" 2>&1
done
