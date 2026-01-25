# ==========================
# RIP-seq Enrichment Analysis (GO, KEGG, Reactome)
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("config.R")

library(tidyverse)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(ReactomePA)

# --------------------------
# Load annotated significant peaks
# --------------------------
WT_sig_annot <- read.csv(file.path(OUTPUT_DIR, "WT_significant_peaks_with_coords.csv"))
A_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "A_significant_peaks_with_coords.csv"))
D_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "D_significant_peaks_with_coords.csv"))

# --------------------------
# Convert Gene IDs to ENTREZ IDs
# --------------------------
convert_to_entrez <- function(gene_ids) {
  entrez <- mapIds(
    org.Sc.sgd.db,
    keys = as.character(gene_ids),
    column = "ENTREZID",
    keytype = GENE_ID_TYPE,
    multiVals = "first"
  )
  unname(entrez[!is.na(entrez)])
}

geneid_WT <- convert_to_entrez(WT_sig_annot$Gene_id)
geneid_A  <- convert_to_entrez(A_sig_annot$Gene_id)
geneid_D  <- convert_to_entrez(D_sig_annot$Gene_id)

# --------------------------
# GO Enrichment
# --------------------------
ego_WT <- enrichGO(
  gene          = geneid_WT,
  OrgDb         = org.Sc.sgd.db,
  ont           = GO_ONTOLOGY,
  pAdjustMethod = PADJ_METHOD,
  qvalueCutoff  = QVAL_CUTOFF
)

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

# --------------------------
# KEGG Enrichment
# --------------------------
ekegg_WT <- enrichKEGG(
  gene          = geneid_WT,
  organism      = KEGG_ORGANISM,
  pAdjustMethod = PADJ_METHOD,
  qvalueCutoff  = QVAL_CUTOFF
)

ekegg_A <- enrichKEGG(
  gene          = geneid_A,
  organism      = KEGG_ORGANISM,
  pAdjustMethod = PADJ_METHOD,
  qvalueCutoff  = QVAL_CUTOFF
)

ekegg_D <- enrichKEGG(
  gene          = geneid_D,
  organism      = KEGG_ORGANISM,
  pAdjustMethod = PADJ_METHOD,
  qvalueCutoff  = QVAL_CUTOFF
)

# --------------------------
# Reactome Enrichment
# --------------------------
ereact_WT <- enrichPathway(
  gene          = geneid_WT,
  organism      = organism,
  pAdjustMethod = PADJ_METHOD,
  qvalueCutoff  = QVAL_CUTOFF
)

ereact_A <- enrichPathway(
  gene          = geneid_A,
  organism      = organism,
  pAdjustMethod = PADJ_METHOD,
  qvalueCutoff  = QVAL_CUTOFF
)

ereact_D <- enrichPathway(
  gene          = geneid_D,
  organism      = organism,
  pAdjustMethod = PADJ_METHOD,
  qvalueCutoff  = QVAL_CUTOFF
)

# --------------------------
# Save results (RDS + CSV)
# --------------------------
saveRDS(ego_WT,   file.path(OUTPUT_DIR, "GO_enrichment_WT.rds"))
saveRDS(ego_A,    file.path(OUTPUT_DIR, "GO_enrichment_A.rds"))
saveRDS(ego_D,    file.path(OUTPUT_DIR, "GO_enrichment_D.rds"))

saveRDS(ekegg_WT, file.path(OUTPUT_DIR, "KEGG_enrichment_WT.rds"))
saveRDS(ekegg_A,  file.path(OUTPUT_DIR, "KEGG_enrichment_A.rds"))
saveRDS(ekegg_D,  file.path(OUTPUT_DIR, "KEGG_enrichment_D.rds"))

saveRDS(ereact_WT, file.path(OUTPUT_DIR, "Reactome_enrichment_WT.rds"))
saveRDS(ereact_A,  file.path(OUTPUT_DIR, "Reactome_enrichment_A.rds"))
saveRDS(ereact_D,  file.path(OUTPUT_DIR, "Reactome_enrichment_D.rds"))

write.csv(as.data.frame(ego_WT),   file.path(OUTPUT_DIR, "GO_enrichment_WT.csv"))
write.csv(as.data.frame(ego_A),    file.path(OUTPUT_DIR, "GO_enrichment_A.csv"))
write.csv(as.data.frame(ego_D),    file.path(OUTPUT_DIR, "GO_enrichment_D.csv"))

write.csv(as.data.frame(ekegg_WT), file.path(OUTPUT_DIR, "KEGG_enrichment_WT.csv"))
write.csv(as.data.frame(ekegg_A),  file.path(OUTPUT_DIR, "KEGG_enrichment_A.csv"))
write.csv(as.data.frame(ekegg_D),  file.path(OUTPUT_DIR, "KEGG_enrichment_D.csv"))

write.csv(as.data.frame(ereact_WT), file.path(OUTPUT_DIR, "Reactome_enrichment_WT.csv"))
write.csv(as.data.frame(ereact_A),  file.path(OUTPUT_DIR, "Reactome_enrichment_A.csv"))
write.csv(as.data.frame(ereact_D),  file.path(OUTPUT_DIR, "Reactome_enrichment_D.csv"))

cat("Enrichment analysis (GO, KEGG, Reactome) completed. Results saved in:", OUTPUT_DIR, "\n")
