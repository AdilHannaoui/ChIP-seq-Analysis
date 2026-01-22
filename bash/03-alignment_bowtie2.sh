#!/bin/bash
set -euo pipefail

# Purpose: Align trimmed RIP-seq reads to the Saccharomyces cerevisiae genome using Bowtie2.
# Output: Sorted and indexed BAM files + BedGraph coverage tracks.
# Author: Adil Hannaoui Anaaoui

WORKDIR="/mnt/c/Users/rna-seq"
BOWTIE2_INDEX="$WORKDIR/BOWTIE2/cerevisiae/index/genome"
OUTPUT_DIR="$WORKDIR/output"
FASTQ_DIR="$OUTPUT_DIR/fastq_trimmed"
THREADS=6

cd "$WORKDIR"

mkdir -p "$OUTPUT_DIR/bowtie2" \
         "$OUTPUT_DIR/visualization"

for FASTQ_FILE in "$FASTQ_DIR"/*.fastq; do
    SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq)
    echo "Processing sample: $SAMPLE_NAME"

    BAM_FILE="$OUTPUT_DIR/bowtie2/${SAMPLE_NAME}_trimmed.bam"
    BEDGRAPH_FILE="$OUTPUT_DIR/visualization/${SAMPLE_NAME}.bedgraph"
    LOG="$OUTPUT_DIR/bowtie2/${SAMPLE_NAME}.log"

    # Run Bowtie2 and sort BAM
    if ! bowtie2 -p "$THREADS" -x "$BOWTIE2_INDEX" -U "$FASTQ_FILE" \
        2> "$LOG" | samtools sort -@ "$THREADS" -o "$BAM_FILE"; then
        echo "Error: Bowtie2 or Samtools failed for $SAMPLE_NAME"
        continue
    fi

    samtools index "$BAM_FILE"

    # Generate BedGraph
    if ! bedtools genomecov -ibam "$BAM_FILE" -bg > "$BEDGRAPH_FILE"; then
        echo "Error: Bedtools failed for $SAMPLE_NAME"
        continue
    fi

    echo "Bowtie2 finished for $SAMPLE_NAME!"
done
