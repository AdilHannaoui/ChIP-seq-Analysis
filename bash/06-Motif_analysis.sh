#!/usr/bin/env bash
# ==========================
# RIP-seq Motif Analysis (MEME Suite)
# Author: Adil Hannaoui Anaaoui
# ==========================

# Exit on error
set -euo pipefail

# --------------------------
# Configuration
# --------------------------
source "$(dirname "$0")/config.sh"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_R="R/config.R"

# Extraer OUTPUT_DIR desde config.R (opcional, si quieres mantener una sola fuente de verdad)
OUTPUT_DIR=$(Rscript -e "source('${CONFIG_R}'); cat(OUTPUT_DIR)")

FASTA_WT="${OUTPUT_DIR}/WT_peak_sequences.fasta"
FASTA_A="${OUTPUT_DIR}/A_peak_sequences.fasta"
FASTA_D="${OUTPUT_DIR}/D_peak_sequences.fasta"

MOTIF_DIR="${OUTPUT_DIR}/motifs"
mkdir -p "${MOTIF_DIR}"

# Parámetros MEME/DREME/FIMO (ajústalos a tu gusto)
NMOTIFS=10
MEME_MINW=6
MEME_MAXW=15
MEME_MOD="zoops"   # o "anr", "oops"

# --------------------------
# Helper function
# --------------------------
run_motif_pipeline () {
    local label="$1"
    local fasta="$2"
    local outdir="${MOTIF_DIR}/${label}"

    echo ">>> Running motif analysis for ${label}"
    mkdir -p "${outdir}"

    # 1) DREME: motivos cortos de novo
    dreme \
        -oc "${outdir}/dreme" \
        -p "${fasta}" \
        -dna \
        -norc

    # 2) MEME: motivos más detallados
    meme "${fasta}" \
        -oc "${outdir}/meme" \
        -dna \
        -mod "${MEME_MOD}" \
        -nmotifs "${NMOTIFS}" \
        -minw "${MEME_MINW}" \
        -maxw "${MEME_MAXW}" \
        -revcomp

    # 3) FIMO: escanear motivos de MEME sobre las mismas secuencias
    fimo \
        --oc "${outdir}/fimo" \
        "${outdir}/meme/meme.xml" \
        "${fasta}"

    echo ">>> Finished motif analysis for ${label}"
    echo
}

# --------------------------
# Run motif analysis per condition
# --------------------------
if [[ ! -f "${FASTA_WT}" ]]; then
    echo "ERROR: FASTA file not found: ${FASTA_WT}" >&2
    exit 1
fi

if [[ ! -f "${FASTA_A}" ]]; then
    echo "ERROR: FASTA file not found: ${FASTA_A}" >&2
    exit 1
fi

if [[ ! -f "${FASTA_D}" ]]; then
    echo "ERROR: FASTA file not found: ${FASTA_D}" >&2
    exit 1
fi

run_motif_pipeline "WT" "${FASTA_WT}"
run_motif_pipeline "A"  "${FASTA_A}"
run_motif_pipeline "D"  "${FASTA_D}"

echo "All motif analyses completed. Results in: ${MOTIF_DIR}"
