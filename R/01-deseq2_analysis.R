# ==========================
# RIP-seq DESeq2 analysis
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("config.R") 

library(DESeq2)
library(tidyverse)

# --------------------------
# Load data
# --------------------------
counts_matrix <- readRDS(COUNTS_MATRIX_PATH)    # counts matrix: genes x samples
colData <- readRDS(SAMPLE_METADATA_PATH)       # metadata: samples x conditions

sample_group <- sub("_.*", "", CONDITIONS)
condition <- sub(".*_", "", CONDITIONS)

colData$sample_group <- factor(sample_group)
colData$condition <- factor(condition, levels = c("IN", "IP"))

colData$group_condition <- factor(paste0(colData$sample_group, "_", colData$condition))

# --------------------------
# Create DESeqDataSet
# --------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = colData,
  design = ~ group_condition
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
# Extract results for contrasts
# --------------------------
res_WT <- results(dds, contrast = c("group_condition", "WT_IP", "WT_IN"))
res_A <- results(dds, contrast = c("group_condition", "Rpb4-S/T-A_IP", "Rpb4-S/T-A_IN"))
res_D <- results(dds, contrast = c("group_condition", "Rpb4-S/T-D_IP", "Rpb4-S/T-D_IN"))

# --------------------------
# Filter significant genes
# --------------------------
res_WT_sig <- res_WT[which(res_WT$padj < PADJ_THRESHOLD), ]
res_A_sig <- res_A[which(res_A$padj < PADJ_THRESHOLD), ]
res_D_sig <- res_D[which(res_D$padj < PADJ_THRESHOLD), ]

# --------------------------
# Save results
# --------------------------
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(dds, file = file.path(OUTPUT_DIR, "dds.rds"))
saveRDS(res_WT_sig, file = file.path(OUTPUT_DIR, "DESeq2_res_WT_sig.rds"))
saveRDS(res_A_sig, file = file.path(OUTPUT_DIR, "DESeq2_res_A_sig.rds"))
saveRDS(res_D_sig, file = file.path(OUTPUT_DIR, "DESeq2_res_D_sig.rds"))

write.csv(as.data.frame(res_WT_sig), file = file.path(OUTPUT_DIR, "DESeq2_res_WT_sig.csv"))
write.csv(as.data.frame(res_A_sig), file = file.path(OUTPUT_DIR, "DESeq2_res_A_sig.csv"))
write.csv(as.data.frame(res_D_sig), file = file.path(OUTPUT_DIR, "DESeq2_res_D_sig.csv"))

cat("DESeq2 analysis completed. Significant results saved in:", OUTPUT_DIR, "\n")
