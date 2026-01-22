#!/bin/bash
set -euo pipefail

# Purpose: Peak calling using MACS2 (robust IP/IN detection)
# Author: Adil Hannaoui Anaaoui

WORKDIR="/mnt/c/Users/rip-seq"
OUTPUT_DIR="$WORKDIR/output"
MAX_REPLICAS=3

cd "$WORKDIR"

mkdir -p "$OUTPUT_DIR/macs2" \
         "$OUTPUT_DIR/logs"

for REP in $(seq 1 "$MAX_REPLICAS"); do

    # Detectar archivos IP con coincidencia insensible a mayúsculas
    for IP_FILE in "$OUTPUT_DIR/bowtie2"/*_[Ii][Pp]${REP}.bam; do
        
        [[ -e "$IP_FILE" ]] || continue  # Si no hay coincidencias, saltar

        BASENAME=$(basename "$IP_FILE")
        
        # Eliminar el sufijo _IP1, _ip1, _Ip1, etc.
        SAMPLE_NAME=$(echo "$BASENAME" | sed -E "s/_[Ii][Pp]${REP}.*//")

        echo "Processing sample: $SAMPLE_NAME (rep $REP)"

        # Buscar archivo IN correspondiente (IN, in, In, iN)
        IN_FILE=$(ls "$OUTPUT_DIR/bowtie2"/"${SAMPLE_NAME}"_[Ii][Nn]${REP}.bam 2>/dev/null || true)

        if [[ -f "$IN_FILE" ]]; then
            
            macs2 callpeak \
                -t "$IP_FILE" \
                -c "$IN_FILE" \
                --format BAM \
                --name "$OUTPUT_DIR/macs2/${SAMPLE_NAME}_rep${REP}" \
                --pvalue 0.1 \
                --call-summits \
                --nomodel \
                --extsize 150 \
                --gsize 1.2e7 \
                > "$OUTPUT_DIR/logs/narrowpeaks_${SAMPLE_NAME}_rep${REP}.log" 2>&1

        else
            echo "ERROR: No IN file found for sample $SAMPLE_NAME (rep $REP)"
        fi
    done
done
