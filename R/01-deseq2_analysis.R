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
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
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

# Extract sample group (M1, M2, M3...) from filenames
sample_groups <- sub("_common_counts\\.txt$", "", basename(COUNTS_FILES))

# Expand sample_groups to match IP/IN columns
expanded_groups <- rep(sample_groups, each = 2)

# Extract condition from column names
conditions <- ifelse(grepl("^IP", sample_names), "IP", "IN")

colData <- data.frame(
  row.names = sample_names,
  sample_group = factor(expanded_groups),
  condition = factor(conditions, levels = c("IN", "IP"))
)

colData$group_condition <- factor(
  paste0(colData$sample_group, "_", colData$condition)
)

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
