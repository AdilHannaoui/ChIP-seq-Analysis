#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq Motif Analysis (MEME Suite, parallel)
# Author: Adil Hannaoui Anaaoui
# ==========================

# Load global config
source "$(dirname "$0")/config.sh"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_R="R/config.R"

OUTPUT_DIR=$(Rscript -e "source('${CONFIG_R}'); cat(OUTPUT_DIR)")

FASTA_WT="${OUTPUT_DIR}/WT_peak_sequences.fasta"
FASTA_A="${OUTPUT_DIR}/A_peak_sequences.fasta"
FASTA_D="${OUTPUT_DIR}/D_peak_sequences.fasta"

MOTIF_DIR="${OUTPUT_DIR}/motifs"
mkdir -p "${MOTIF_DIR}"
mkdir -p "${OUTPUT_DIR}/logs"

# MEME/DREME/FIMO parameters
NMOTIFS=10
MEME_MINW=6
MEME_MAXW=15
MEME_MOD="zoops"

#############################################
# Function to run motif analysis for one label
#############################################
run_motif_pipeline() {
    local label="$1"
    local fasta="$2"

    local outdir="${MOTIF_DIR}/${label}"
    local logfile="${OUTPUT_DIR}/logs/motif_${label}.log"

    echo ">>> Running motif analysis for ${label}"

    {
        mkdir -p "${outdir}"

        # 1) DREME
        dreme \
            -oc "${outdir}/dreme" \
            -p "${fasta}" \
            -dna \
            -norc

        # 2) MEME
        meme "${fasta}" \
            -oc "${outdir}/meme" \
            -dna \
            -mod "${MEME_MOD}" \
            -nmotifs "${NMOTIFS}" \
            -minw "${MEME_MINW}" \
            -maxw "${MEME_MAXW}" \
            -revcomp

        # 3) FIMO
        fimo \
            --oc "${outdir}/fimo" \
            "${outdir}/meme/meme.xml" \
            "${fasta}"

        echo ">>> Finished motif analysis for ${label}"
        echo

    } > "${logfile}" 2>&1
}

export -f run_motif_pipeline
export MOTIF_DIR OUTPUT_DIR NMOTIFS MEME_MINW MEME_MAXW MEME_MOD

#############################################
# Validate FASTA files
#############################################
for f in "$FASTA_WT" "$FASTA_A" "$FASTA_D"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: FASTA file not found: $f" >&2
        exit 1
    fi
done

#############################################
# Run motif analysis in parallel
#############################################
parallel -j "$THREADS" run_motif_pipeline ::: \
    "WT" "A" "D" ::: \
    "$FASTA_WT" "$FASTA_A" "$FASTA_D"

echo "All motif analyses completed. Results in: ${MOTIF_DIR}"
