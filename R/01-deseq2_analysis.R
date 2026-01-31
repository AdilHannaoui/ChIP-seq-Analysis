# ==========================
# RIP-seq DESeq2 analysis
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

library(DESeq2)
library(dplyr)

# --------------------------
# Load data
# --------------------------
COUNTS_FILES <- list.files(
  path = "output/macs2",
  pattern = "_common_counts\\.txt$",
  full.names = TRUE
)

count_list <- lapply(COUNTS_FILES, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_name <- sub("_common_counts\\.txt$", "", basename(f))  # M1, M12, WT
  colnames(df)[9:14] <- paste0(sample_name, "_", colnames(df)[9:14])
  df
})

# Merge all count tables by genomic annotation columns
counts_merged <- Reduce(function(x, y) merge(
  x, y,
  by = c("chrom","start","end","name","score","strand","Gene_id","Biotype"),
  all = TRUE
), count_list)

# --------------------------
# Collapse peaks by gene (SUM)
# --------------------------
count_data <- counts_merged %>%
  select(-(chrom:Biotype)) %>%
  mutate(Gene_id = counts_merged$Gene_id) %>%
  group_by(Gene_id) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

counts_matrix <- as.data.frame(count_data)
rownames(counts_matrix) <- counts_matrix$Gene_id
counts_matrix$Gene_id <- NULL

# --------------------------
# Build colData automatically
# --------------------------
sample_names <- colnames(counts_matrix)

sample_group <- sub("_.*", "", sample_names)                 # M1, M12, WT
condition    <- ifelse(grepl("_IP", sample_names), "IP", "IN")

colData <- data.frame(
  row.names     = sample_names,
  sample_group  = factor(sample_group),
  condition     = factor(condition, levels = c("IN", "IP"))
)

# --------------------------
# Create DESeqDataSet
# --------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = colData,
  design    = ~ sample_group + condition + sample_group:condition
)

# --------------------------
# Filter low-count genes
# --------------------------
dds <- dds[rowSums(counts(dds)) > MIN_COUNTS_FILTER, ]

# --------------------------
# Run DESeq2
# --------------------------
dds <- DESeq(dds)

# --------------------------
# Extract results: IP vs IN per sample
# --------------------------
res_M1 <- results(dds, name = "condition_IP_vs_IN")
res_M12 <- results(dds, name = "sample_groupM12.conditionIP")
res_WT  <- results(dds, name = "sample_groupWT.conditionIP")

res_M1_sig  <- res_M1[which(res_M1$padj  < PADJ_THRESHOLD), ]
res_M12_sig <- res_M12[which(res_M12$padj < PADJ_THRESHOLD), ]
res_WT_sig  <- res_WT[which(res_WT$padj  < PADJ_THRESHOLD), ]

# --------------------------
# Extract results: differences between groups
# --------------------------
res_WT_vs_M1 <- results(dds, contrast = list(
  c("sample_groupWT.conditionIP"),
  c("sample_groupM1.conditionIP")
))

res_WT_vs_M12 <- results(dds, contrast = list(
  c("sample_groupWT.conditionIP"),
  c("sample_groupM12.conditionIP")
))

res_M1_vs_M12 <- results(dds, contrast = list(
  c("sample_groupM1.conditionIP"),
  c("sample_groupM12.conditionIP")
))

res_WT_vs_M1_sig  <- res_WT_vs_M1[which(res_WT_vs_M1$padj  < PADJ_THRESHOLD), ]
res_WT_vs_M12_sig <- res_WT_vs_M12[which(res_WT_vs_M12$padj < PADJ_THRESHOLD), ]
res_M1_vs_M12_sig <- res_M1_vs_M12[which(res_M1_vs_M12$padj < PADJ_THRESHOLD), ]

# --------------------------
# Save results
# --------------------------
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Core objects for downstream modules
saveRDS(dds,           file = file.path(OUTPUT_DIR, "dds.rds"))
saveRDS(colData,       file = file.path(OUTPUT_DIR, "colData.rds"))
saveRDS(counts_merged, file = file.path(OUTPUT_DIR, "counts_merged.rds"))

# IP vs IN per sample
saveRDS(res_M1_sig,  file = file.path(OUTPUT_DIR, "DESeq2_res_M1_sig.rds"))
saveRDS(res_M12_sig, file = file.path(OUTPUT_DIR, "DESeq2_res_M12_sig.rds"))
saveRDS(res_WT_sig,  file = file.path(OUTPUT_DIR, "DESeq2_res_WT_sig.rds"))

write.csv(as.data.frame(res_M1_sig),  file = file.path(OUTPUT_DIR, "DESeq2_res_M1_sig.csv"))
write.csv(as.data.frame(res_M12_sig), file = file.path(OUTPUT_DIR, "DESeq2_res_M12_sig.csv"))
write.csv(as.data.frame(res_WT_sig),  file = file.path(OUTPUT_DIR, "DESeq2_res_WT_sig.csv"))

# Differences between groups
saveRDS(res_WT_vs_M1_sig,  file = file.path(OUTPUT_DIR, "DESeq2_res_WT_vs_M1_sig.rds"))
saveRDS(res_WT_vs_M12_sig, file = file.path(OUTPUT_DIR, "DESeq2_res_WT_vs_M12_sig.rds"))
saveRDS(res_M1_vs_M12_sig, file = file.path(OUTPUT_DIR, "DESeq2_res_M1_vs_M12_sig.rds"))

write.csv(as.data.frame(res_WT_vs_M1_sig),  file = file.path(OUTPUT_DIR, "DESeq2_res_WT_vs_M1_sig.csv"))
write.csv(as.data.frame(res_WT_vs_M12_sig), file = file.path(OUTPUT_DIR, "DESeq2_res_WT_vs_M12_sig.csv"))
write.csv(as.data.frame(res_M1_vs_M12_sig), file = file.path(OUTPUT_DIR, "DESeq2_res_M1_vs_M12_sig.csv"))

cat("DESeq2 analysis completed. Significant results saved in:", OUTPUT_DIR, "\n")
