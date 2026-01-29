# ==========================
# RIP-seq Peak Coordinates Extraction
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

library(tidyverse)

# --------------------------
# Load DESeq2 significant results
# --------------------------
res_WT_sig <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_WT_sig.rds"))
res_A_sig  <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_A_sig.rds"))
res_D_sig  <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_D_sig.rds"))

# --------------------------
# Load counts with coordinates
# --------------------------
counts_df <- read.delim(COUNTS_MATRIX_PATH, header = TRUE)

# Ensure peak ID is used as rownames
rownames(counts_df) <- counts_df$name

# --------------------------
# Merge coordinates with DESeq2 results
# --------------------------
WT_sig_annot <- cbind(
  counts_df[rownames(res_WT_sig), c("chrom","start","end")],
  as.data.frame(res_WT_sig)
)

A_sig_annot <- cbind(
  counts_df[rownames(res_A_sig), c("chrom","start","end")],
  as.data.frame(res_A_sig)
)

D_sig_annot <- cbind(
  counts_df[rownames(res_D_sig), c("chrom","start","end")],
  as.data.frame(res_D_sig)
)

# --------------------------
# Save annotated CSV tables
# --------------------------
write.csv(WT_sig_annot, file.path(OUTPUT_DIR, "WT_significant_peaks_with_coords.csv"), row.names = TRUE)
write.csv(A_sig_annot,  file.path(OUTPUT_DIR, "A_significant_peaks_with_coords.csv"),  row.names = TRUE)
write.csv(D_sig_annot,  file.path(OUTPUT_DIR, "D_significant_peaks_with_coords.csv"),  row.names = TRUE)

# --------------------------
# Export BED files for IGV
# --------------------------
WT_bed <- WT_sig_annot %>% select(chrom, start, end)
A_bed  <- A_sig_annot  %>% select(chrom, start, end)
D_bed  <- D_sig_annot  %>% select(chrom, start, end)

write.table(WT_bed, file.path(OUTPUT_DIR, "WT_significant_peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(A_bed,  file.path(OUTPUT_DIR, "A_significant_peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(D_bed,  file.path(OUTPUT_DIR, "D_significant_peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Peak coordinate extraction completed. CSV and BED files saved in:", OUTPUT_DIR, "\n")

