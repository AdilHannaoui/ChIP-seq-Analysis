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
FASTA_M1="${OUTPUT_DIR}/M1_peak_sequences.fasta"
FASTA_M12="${OUTPUT_DIR}/M12_peak_sequences.fasta"
FASTA_M1_M12="${OUTPUT_DIR}/M1_M12_peak_sequences.fasta"
FASTA_M12_WT="${OUTPUT_DIR}/M12_WT_peak_sequences.fasta"

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

    echo ">>> Running motif analysis for ${label} (fasta: ${fasta})"

    {
        mkdir -p "${outdir}"

        # Create a temporary FASTA with unique headers for MEME (preserve sequences)
        tmp_fasta="${outdir}/input_for_meme.fasta"
        awk -v lbl="${label}" 'BEGIN{n=0}
            /^>/{ n++; print ">" lbl "_" n; next }
            { print }
        ' "${fasta}" > "${tmp_fasta}"

        # Count sequences in the tmp fasta
        seq_count=$(grep -c "^>" "${tmp_fasta}" || true)
        echo "INFO: sequences in ${tmp_fasta}: ${seq_count}"

        # If fewer than 2 sequences, skip MEME (unless you want -mod anr)
        if [[ "${seq_count}" -lt 2 ]]; then
            echo "WARNING: Not enough sequences for MEME (need >=2). Skipping MEME and FIMO for ${label}."
            # Still run DREME if you want (DREME can work with single-sequence cases differently)
            dreme -oc "${outdir}/dreme" -p "${fasta}" -dna -norc || echo "DREME failed or produced no output"
            echo ">>> Finished motif analysis for ${label} (skipped MEME/FIMO)"
            echo
            continue 2  # exit the subshell block and finish
        fi

        # 1) DREME (run on original fasta or tmp_fasta as you prefer)
        dreme -oc "${outdir}/dreme" -p "${fasta}" -dna -norc

        # 2) MEME (run on tmp_fasta with unique headers)
        meme "${tmp_fasta}" \
            -oc "${outdir}/meme" \
            -dna \
            -mod "${MEME_MOD}" \
            -nmotifs "${NMOTIFS}" \
            -minw "${MEME_MINW}" \
            -maxw "${MEME_MAXW}" \
            -revcomp \
            -maxsize 1000000 -p 1

        # 3) FIMO (only if meme.xml exists and non-empty)
        if [[ -s "${outdir}/meme/meme.xml" ]]; then
            fimo --oc "${outdir}/fimo" "${outdir}/meme/meme.xml" "${fasta}"
        else
            echo "No meme.xml or empty for ${label}, skipping FIMO"
        fi

        echo ">>> Finished motif analysis for ${label}"
        echo

    } > "${logfile}" 2>&1
}

export -f run_motif_pipeline
export MOTIF_DIR OUTPUT_DIR NMOTIFS MEME_MINW MEME_MAXW MEME_MOD

#############################################
# Validate FASTA files and build pairs
# Headers may contain a dot (.) — allow any header content
#############################################
labels=(WT M1 M12 M1_M12 M12_WT)
fastas=("$FASTA_WT" "$FASTA_M1" "$FASTA_M12" "$FASTA_M1_M12" "$FASTA_M12_WT")

# Check files exist and have sequences; build arrays of valid pairs
valid_labels=()
valid_fastas=()

for i in "${!labels[@]}"; do
    lbl="${labels[$i]}"
    f="${fastas[$i]}"

    if [[ ! -f "$f" ]]; then
        echo "WARNING: FASTA not found for ${lbl}: ${f} (skipping)"
        continue
    fi

    # Count sequences (headers starting with >)
    seq_count=$(grep -c "^>" "$f" || true)
    if [[ "$seq_count" -eq 0 ]]; then
        echo "WARNING: FASTA for ${lbl} has 0 sequences: ${f} (skipping)"
        continue
    fi

    # Validate sequence characters only (ignore header lines).
    # Support multi-line sequences: join sequence lines and check characters ACGTN only.
    # We'll scan file and print first offending line if any.
    bad_seq=$(awk 'BEGIN{bad=0}
        /^>/{next}
        {
          if($0 !~ /^[ACGTNacgtn]+$/) {
            print FILENAME ":" NR ":" $0
            bad=1
            exit
          }
        }
        END{ if(bad==1) exit 1 }' "$f" 2>/dev/null || true)

    if [[ -n "$bad_seq" ]]; then
        echo "WARNING: FASTA ${f} contains unexpected characters in sequence lines (first bad line shown):"
        echo "$bad_seq"
        echo "Skipping ${lbl}"
        continue
    fi

    echo "OK: ${lbl} -> ${f} (sequences: ${seq_count})"
    valid_labels+=("$lbl")
    valid_fastas+=("$f")
done

if [[ "${#valid_labels[@]}" -eq 0 ]]; then
    echo "ERROR: No valid FASTA files found. Exiting."
    exit 1
fi

#############################################
# Run motif analysis in parallel (1:1 pairing)
#############################################
# Use --link to ensure 1:1 pairing between labels and fastas
# THREADS should be defined in config.sh; fallback to 1 if not set
THREADS="${THREADS:-1}"
parallel -j "${THREADS}" --link run_motif_pipeline ::: "${valid_labels[@]}" ::: "${valid_fastas[@]}"

echo "All motif analyses completed. Results in: ${MOTIF_DIR}"
