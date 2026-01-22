#!/bin/bash
set -euo pipefail

# Purpose: Obtain intersection peaks across all replicates using bedtools multiinter
# Author: Adil Hannaoui Anaaoui

WORKDIR="/mnt/c/Users/rip-seq"
OUTPUT_DIR="$WORKDIR/output"
GTF_FILE="$WORKDIR/annotation/Saccharomyces_cerevisiae.R64-1-1.112.gtf"

mkdir -p "$OUTPUT_DIR/macs2" "$OUTPUT_DIR/logs"

#############################################
# 1. Detect Samples
#############################################

for FILE in "$OUTPUT_DIR/macs2"/*_rep1_peaks.narrowPeak; do
    [[ -e "$FILE" ]] || continue

    BASENAME=$(basename "$FILE" | sed 's/_rep1_peaks.narrowPeak//')
    echo "Processing sample: $BASENAME"

    #############################################
    # 2. Collect Peak Files
    #############################################

    NARROWPEAK_FILES=()
    for NP in "$OUTPUT_DIR/macs2/${BASENAME}"_rep*_peaks.narrowPeak; do
        [[ -e "$NP" ]] || continue
        NARROWPEAK_FILES+=("$NP")
    done

    NUM_REPS=${#NARROWPEAK_FILES[@]}

    if [[ $NUM_REPS -lt 2 ]]; then
        echo "Not enough replicates for $BASENAME"
        continue
    fi

    echo "Found $NUM_REPS replicates for $BASENAME"

    #############################################
    # 3. (Common peaks), Multiple Intersection
    #############################################

    bedtools multiinter -i "${NARROWPEAK_FILES[@]}" \
        | awk -v n="$NUM_REPS" '$4 == n {print $1, $2, $3}' OFS="\t" \
        > "$OUTPUT_DIR/macs2/${BASENAME}_common_raw.bed"

    #############################################
    # 4. GTF Annotation
    #############################################

    annotate() {
        bedtools intersect -a "$1" -b "$GTF_FILE" -wa -wb |
        awk -F'\t' '
        {
            gene_id="NA"; biotype="NA"
            if (match($0, /gene_id "([^"]+)"/, m)) gene_id=m[1]
            if (match($0, /gene_biotype "([^"]+)"/, b)) biotype=b[1]
            print $1, $2, $3, ".", ".", ".", gene_id, biotype
        }' OFS="\t" | sort -u
    }

    annotate "$OUTPUT_DIR/macs2/${BASENAME}_common_raw.bed" \
        > "$OUTPUT_DIR/macs2/${BASENAME}_common_annotated.bed"

    rm "$OUTPUT_DIR/macs2/${BASENAME}_common_raw.bed"

    #############################################
    # 5. Header Obtaining
    #############################################

    HEADER="chrom\tstart\tend\tname\tscore\tstrand\tGene_id\tBiotype"

    for REP in $(seq 1 "$NUM_REPS"); do
        HEADER+="\tIP${REP}\tIN${REP}"
    done

    #############################################
    # 6. Collect BAMs
    #############################################

    BAM_FILES=()
    for REP in $(seq 1 "$NUM_REPS"); do
        BAM_FILES+=("$OUTPUT_DIR/bowtie2/${BASENAME}_IP${REP}.bam")
        BAM_FILES+=("$OUTPUT_DIR/bowtie2/${BASENAME}_IN${REP}.bam")
    done

    #############################################
    # 7. multicov
    #############################################

    bedtools multicov -bams "${BAM_FILES[@]}" \
        -bed "$OUTPUT_DIR/macs2/${BASENAME}_common_annotated.bed" \
        > "$OUTPUT_DIR/macs2/${BASENAME}_common_counts.tmp"

    #############################################
    # 8. Combine header + counts
    #############################################

    {
        echo -e "$HEADER"
        cat "$OUTPUT_DIR/macs2/${BASENAME}_common_counts.tmp"
    } > "$OUTPUT_DIR/macs2/${BASENAME}_common_counts.txt"

    rm "$OUTPUT_DIR/macs2/${BASENAME}_common_counts.tmp"

    echo "Finished sample: $BASENAME"
done
