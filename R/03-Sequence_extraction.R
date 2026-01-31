# ==========================
# RIP-seq Significant Peak Sequence Extraction
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

genome <- REFERENCE_GENOME

# --------------------------
# Helper: normalize chromosome names
# --------------------------
normalize_chr <- function(x) {
  x <- as.character(x)
  x[x == "Mito"] <- "M"        # Mito → M
  paste0("chr", x)             # agregar prefijo chr
}

# --------------------------
# Helper: safe GRanges + FASTA extraction
# --------------------------
process_peaks <- function(df, label, fasta_name) {
  if (is.null(df) || nrow(df) == 0) {
    message("No significant peaks for ", label, " (0 rows). Skipping FASTA generation.")
    return(NULL)
  }
  if (!all(c("chrom", "start", "end", "name") %in% colnames(df))) {
    stop(paste0(
      "Data frame for ", label,
      " does not have required columns: chrom, start, end, name. Got: ",
      paste(colnames(df), collapse = ", ")
    ))
  }
  
  df$chrom <- normalize_chr(df$chrom)
  
  gr <- GRanges(
    seqnames = df$chrom,
    ranges   = IRanges(start = df$start, end = df$end)
  )
  
  seqs <- getSeq(genome, gr)
  names(seqs) <- df$name
  
  out_path <- file.path(OUTPUT_DIR, fasta_name)
  writeXStringSet(seqs, out_path)
  message("FASTA written for ", label, " → ", out_path)
  
  invisible(seqs)
}

# --------------------------
# Load annotated significant peaks (CSV)
# --------------------------
WT_sig_annot      <- read.csv(file.path(OUTPUT_DIR, "WT_significant_peaks_with_coords.csv"))
M1_sig_annot      <- read.csv(file.path(OUTPUT_DIR, "M1_significant_peaks_with_coords.csv"))
M12_sig_annot     <- read.csv(file.path(OUTPUT_DIR, "M12_significant_peaks_with_coords.csv"))

WT_M1_sig_annot   <- read.csv(file.path(OUTPUT_DIR, "WT_M1_significant_peaks_with_coords.csv"))
M1_M12_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "M1_M12_significant_peaks_with_coords.csv"))
M12_WT_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "M12_WT_significant_peaks_with_coords.csv"))

# --------------------------
# Process each contrast
# --------------------------
seqs_WT      <- process_peaks(WT_sig_annot,     "WT",      "WT_peak_sequences.fasta")
seqs_M1      <- process_peaks(M1_sig_annot,     "M1",      "M1_peak_sequences.fasta")
seqs_M12     <- process_peaks(M12_sig_annot,    "M12",     "M12_peak_sequences.fasta")

seqs_WT_M1   <- process_peaks(WT_M1_sig_annot,  "WT_M1",   "WT_M1_peak_sequences.fasta")
seqs_M1_M12  <- process_peaks(M1_M12_sig_annot, "M1_M12",  "M1_M12_peak_sequences.fasta")
seqs_M12_WT  <- process_peaks(M12_WT_sig_annot, "M12_WT",  "M12_WT_peak_sequences.fasta")

cat("FASTA extraction module finished. Check messages above for skipped contrasts.\n")

