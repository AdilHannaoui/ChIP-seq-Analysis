# ==========================
# RIP-seq Peak Coordinates Extraction
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

library(dplyr)

# --------------------------
# Load DESeq2 significant results
# --------------------------
res_WT_sig <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_WT_sig.rds"))
res_M1_sig  <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_M1_sig.rds"))
res_M12_sig  <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_M12_sig.rds"))

res_WT_vs_M1_sig <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_WT_vs_M1_sig.rds"))
res_M12_vs_M1_sig  <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_M12_vs_M1_sig.rds"))
res_WT_vs_M12_sig  <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_WT_vs_M12_sig.rds"))

# --------------------------
# Load counts with coordinates
# --------------------------
counts_df <- readRDS(file.path(OUTPUT_DIR, "counts_merged.rds"))
colnames(counts_df)


# --------------------------
# Merge individuals coordinates with DESeq2 results
# --------------------------
WT_sig_annot <- counts_df %>%
  filter(Gene_id %in% rownames(res_WT_sig)) %>%
  left_join(
    as.data.frame(res_WT_sig) %>% mutate(Gene_id = rownames(res_WT_sig)),
    by = "Gene_id"
  )

M1_sig_annot <- counts_df %>%
  filter(Gene_id %in% rownames(res_M1_sig)) %>%
  left_join(
    as.data.frame(res_M1_sig) %>% mutate(Gene_id = rownames(res_M1_sig)),
    by = "Gene_id"
  )

M12_sig_annot <- counts_df %>%
  filter(Gene_id %in% rownames(res_M12_sig)) %>%
  left_join(
    as.data.frame(res_M12_sig) %>% mutate(Gene_id = rownames(res_M12_sig)),
    by = "Gene_id"
  )

# --------------------------
# Merge paired coordinates with DESeq2 results
# --------------------------
WT_M1_sig_annot <- counts_df %>%
  filter(Gene_id %in% rownames(res_WT_vs_M1_sig)) %>%
  left_join(
    as.data.frame(res_WT_vs_M1_sig) %>% mutate(Gene_id = rownames(res_WT_vs_M1_sig)),
    by = "Gene_id"
  )

M1_M12_sig_annot <- counts_df %>%
  filter(Gene_id %in% rownames(res_M12_vs_M1_sig)) %>%
  left_join(
    as.data.frame(res_M12_vs_M1_sig) %>% mutate(Gene_id = rownames(res_M12_vs_M1_sig)),
    by = "Gene_id"
  )

M12_WT_sig_annot <- counts_df %>%
  filter(Gene_id %in% rownames(res_WT_vs_M12_sig)) %>%
  left_join(
    as.data.frame(res_WT_vs_M12_sig) %>% mutate(Gene_id = rownames(res_WT_vs_M12_sig)),
    by = "Gene_id"
  )

# --------------------------
# Save annotated CSV tables
# --------------------------
write.csv(WT_sig_annot, file.path(OUTPUT_DIR, "WT_significant_peaks_with_coords.csv"), row.names = TRUE)
write.csv(M1_sig_annot,  file.path(OUTPUT_DIR, "M1_significant_peaks_with_coords.csv"),  row.names = TRUE)
write.csv(M12_sig_annot,  file.path(OUTPUT_DIR, "M12_significant_peaks_with_coords.csv"),  row.names = TRUE)

write.csv(WT_M1_sig_annot, file.path(OUTPUT_DIR, "WT_M1_significant_peaks_with_coords.csv"), row.names = TRUE)
write.csv(M1_M12_sig_annot,  file.path(OUTPUT_DIR, "M1_M12_significant_peaks_with_coords.csv"),  row.names = TRUE)
write.csv(M12_WT_sig_annot,  file.path(OUTPUT_DIR, "M12_WT_significant_peaks_with_coords.csv"),  row.names = TRUE)
# --------------------------
# Export BED files for IGV
# --------------------------
WT_bed <- WT_sig_annot %>% select(chrom, start, end)
M1_bed  <- M1_sig_annot  %>% select(chrom, start, end)
M12_bed  <- M12_sig_annot  %>% select(chrom, start, end)

WT_M1_bed <- WT_M1_sig_annot %>% select(chrom, start, end)
M1_M12_bed  <- M1_M12_sig_annot  %>% select(chrom, start, end)
M12_WT_bed  <- M12_WT_sig_annot  %>% select(chrom, start, end)


write.table(WT_bed, file.path(OUTPUT_DIR, "WT_significant_peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(M1_bed,  file.path(OUTPUT_DIR, "M1_significant_peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(M12_bed,  file.path(OUTPUT_DIR, "M12_significant_peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



write.table(WT_M1_bed, file.path(OUTPUT_DIR, "WT_M1_significant_peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(M1_M12_bed,  file.path(OUTPUT_DIR, "M1_M12_significant_peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(M12_WT_bed,  file.path(OUTPUT_DIR, "M12_WT_significant_peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Peak coordinate extraction completed. CSV and BED files saved in:", OUTPUT_DIR, "\n")
