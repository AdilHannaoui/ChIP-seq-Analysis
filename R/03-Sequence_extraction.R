# ==========================
# RIP-seq Significant Peak Sequence Extraction
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

library(tidyverse)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(BSgenome.Scerevisiae.UCSC.sacCer3)  

genome <- REFERENCE_GENOME

# --------------------------
# Load annotated significant peaks (CSV)
# --------------------------
WT_sig_annot <- read.csv(file.path(OUTPUT_DIR, "WT_significant_peaks_with_coords.csv"))
A_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "A_significant_peaks_with_coords.csv"))
D_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "D_significant_peaks_with_coords.csv"))

# --------------------------
# Build GRanges objects
# --------------------------
peaks_WT <- GRanges(
  seqnames = WT_sig_annot$chrom,
  ranges   = IRanges(start = WT_sig_annot$start, end = WT_sig_annot$end)
)

peaks_A <- GRanges(
  seqnames = A_sig_annot$chrom,
  ranges   = IRanges(start = A_sig_annot$start, end = A_sig_annot$end)
)

peaks_D <- GRanges(
  seqnames = D_sig_annot$chrom,
  ranges   = IRanges(start = D_sig_annot$start, end = D_sig_annot$end)
)

# --------------------------
# Extract FASTA sequences
# --------------------------
seqs_WT <- getSeq(genome, peaks_WT)
seqs_A  <- getSeq(genome, peaks_A)
seqs_D  <- getSeq(genome, peaks_D)

# Add peak names to FASTA headers
names(seqs_WT) <- WT_sig_annot$name
names(seqs_A)  <- A_sig_annot$name
names(seqs_D)  <- D_sig_annot$name

# --------------------------
# Save FASTA files
# --------------------------
writeXStringSet(seqs_WT, file.path(OUTPUT_DIR, "WT_peak_sequences.fasta"))
writeXStringSet(seqs_A,  file.path(OUTPUT_DIR, "A_peak_sequences.fasta"))
writeXStringSet(seqs_D,  file.path(OUTPUT_DIR, "D_peak_sequences.fasta"))

cat("FASTA extraction completed. Files saved in:", OUTPUT_DIR, "\n")
