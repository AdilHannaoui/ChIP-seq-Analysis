# ==========================
# RIP-seq Enrichment Peak Analysis
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("config.R")

library(tidyverse)
library(clusterProfiler)
library(org.Sc.sgd.db)

# --------------------------
# Load annotated significant peaks
# --------------------------
WT_sig_annot <- read.csv(file.path(OUTPUT_DIR, "WT_significant_peaks_with_coords.csv"))
A_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "A_significant_peaks_with_coords.csv"))
D_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "D_significant_peaks_with_coords.csv"))

# --------------------------
# Map gene IDs to ENTREZ
# --------------------------
geneid_A <- mapIds(
  org.Sc.sgd.db,
  keys = as.character(A_sig_annot$Gene_id),
  column = "ENTREZID",
  keytype = GENE_ID_TYPE,
  multiVals = "first"
)

geneid_D <- mapIds(
  org.Sc.sgd.db,
  keys = as.character(D_sig_annot$Gene_id),
  column = "ENTREZID",
  keytype = GENE_ID_TYPE,
  multiVals = "first"
)

geneid_WT <- mapIds(
  org.Sc.sgd.db,
  keys = as.character(WT_sig_annot$Gene_id),
  column = "ENTREZID",
  keytype = GENE_ID_TYPE,
  multiVals = "first"
)

# Remove NAs
geneid_A  <- unname(geneid_A[!is.na(geneid_A)])
geneid_D  <- unname(geneid_D[!is.na(geneid_D)])
geneid_WT <- unname(geneid_WT[!is.na(geneid_WT)])

# --------------------------
# GO Enrichment
# --------------------------
ego_A <- enrichGO(
  gene          = geneid_A,
  OrgDb         = org.Sc.sgd.db,
  ont           = GO_ONTOLOGY,
  pAdjustMethod = PADJ_METHOD,
  qvalueCutoff  = QVAL_CUTOFF
)

ego_D <- enrichGO(
  gene          = geneid_D,
  OrgDb         = org.Sc.sgd.db,
  ont           = GO_ONTOLOGY,
  pAdjustMethod = PADJ_METHOD,
  qvalueCutoff  = QVAL_CUTOFF
)

ego_WT <- enrichGO(
  gene          = geneid_WT,
  OrgDb         = org.Sc.sgd.db,
  ont           = GO_ONTOLOGY,
  pAdjustMethod = PADJ_METHOD,
  qvalueCutoff  = QVAL_CUTOFF
)

# --------------------------
# Save results
# --------------------------
saveRDS(ego_A,  file.path(OUTPUT_DIR, "GO_enrichment_A.rds"))
saveRDS(ego_D,  file.path(OUTPUT_DIR, "GO_enrichment_D.rds"))
saveRDS(ego_WT, file.path(OUTPUT_DIR, "GO_enrichment_WT.rds"))

write.csv(as.data.frame(ego_A),  file.path(OUTPUT_DIR, "GO_enrichment_A.csv"))
write.csv(as.data.frame(ego_D),  file.path(OUTPUT_DIR, "GO_enrichment_D.csv"))
write.csv(as.data.frame(ego_WT), file.path(OUTPUT_DIR, "GO_enrichment_WT.csv"))

cat("GO enrichment analysis completed. Results saved in:", OUTPUT_DIR, "\n")
