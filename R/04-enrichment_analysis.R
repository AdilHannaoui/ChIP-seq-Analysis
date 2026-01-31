# ==========================
# RIP-seq Enrichment Analysis (GO, KEGG, Reactome)
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

library(clusterProfiler)
library(org.Sc.sgd.db)
library(ReactomePA)

# --------------------------
# Helper: safe Entrez conversion
# --------------------------
convert_to_entrez <- function(gene_ids, label) {
  if (is.null(gene_ids) || length(gene_ids) == 0) {
    message("No gene IDs for ", label, ". Skipping Entrez conversion.")
    return(character(0))
  }
  
  entrez <- mapIds(
    org.Sc.sgd.db,
    keys = as.character(gene_ids),
    column = "ENTREZID",
    keytype = GENE_ID_TYPE,
    multiVals = "first"
  )
  
  entrez <- unname(entrez[!is.na(entrez)])
  
  if (length(entrez) == 0) {
    message("No valid Entrez IDs for ", label, ".")
  }
  
  return(entrez)
}

# --------------------------
# Helper: safe enrichment wrapper
# --------------------------
run_enrichment <- function(entrez_ids, orf_ids, label) {
  
  # GO + Reactome need ENTREZ
  if (length(entrez_ids) == 0) {
    message("No Entrez IDs for ", label, ". Skipping GO and Reactome.")
    ego <- NULL
    ereact <- NULL
  } else {
    message("Running GO + Reactome for ", label, " (", length(entrez_ids), " genes)...")
    
    ego <- enrichGO(
      gene          = entrez_ids,
      OrgDb         = org.Sc.sgd.db,
      ont           = GO_ONTOLOGY,
      pAdjustMethod = PADJ_METHOD,
      qvalueCutoff  = QVAL_CUTOFF
    )
    
    ereact <- enrichPathway(
      gene          = entrez_ids,
      organism      = organism,
      pAdjustMethod = PADJ_METHOD,
      qvalueCutoff  = QVAL_CUTOFF
    )
  }
  
  # KEGG needs ORF IDs
  if (length(orf_ids) == 0) {
    message("No ORF IDs for ", label, ". Skipping KEGG.")
    ekegg <- NULL
  } else {
    message("Running KEGG for ", label, " (", length(orf_ids), " ORFs)...")
    
    ekegg <- enrichKEGG(
      gene          = orf_ids,
      organism      = "sce",
      keyType       = "kegg",
      pAdjustMethod = PADJ_METHOD,
      qvalueCutoff  = QVAL_CUTOFF
    )
  }
  
  return(list(GO=ego, KEGG=ekegg, REACT=ereact))
}

# --------------------------
# Load annotated significant peaks
# --------------------------
load_df <- function(path, label) {
  df <- read.csv(path)
  if (nrow(df) == 0) {
    message("No significant peaks for ", label, ". Skipping.")
  }
  return(df)
}

WT_sig_annot     <- load_df(file.path(OUTPUT_DIR, "WT_significant_peaks_with_coords.csv"), "WT")
M1_sig_annot     <- load_df(file.path(OUTPUT_DIR, "M1_significant_peaks_with_coords.csv"), "M1")
M12_sig_annot    <- load_df(file.path(OUTPUT_DIR, "M12_significant_peaks_with_coords.csv"), "M12")

WT_M1_sig_annot  <- load_df(file.path(OUTPUT_DIR, "WT_M1_significant_peaks_with_coords.csv"), "WT_M1")
M1_M12_sig_annot <- load_df(file.path(OUTPUT_DIR, "M1_M12_significant_peaks_with_coords.csv"), "M1_M12")
M12_WT_sig_annot <- load_df(file.path(OUTPUT_DIR, "M12_WT_significant_peaks_with_coords.csv"), "M12_WT")

# --------------------------
# Convert to Entrez (GO + Reactome)
# --------------------------
geneid_WT      <- convert_to_entrez(WT_sig_annot$Gene_id, "WT")
geneid_M1      <- convert_to_entrez(M1_sig_annot$Gene_id, "M1")
geneid_M12     <- convert_to_entrez(M12_sig_annot$Gene_id, "M12")

geneid_WT_M1   <- convert_to_entrez(WT_M1_sig_annot$Gene_id, "WT_M1")
geneid_M1_M12  <- convert_to_entrez(M1_M12_sig_annot$Gene_id, "M1_M12")
geneid_M12_WT  <- convert_to_entrez(M12_WT_sig_annot$Gene_id, "M12_WT")

# --------------------------
# ORF IDs for KEGG
# --------------------------
orf_WT      <- WT_sig_annot$Gene_id
orf_M1      <- M1_sig_annot$Gene_id
orf_M12     <- M12_sig_annot$Gene_id

orf_WT_M1   <- WT_M1_sig_annot$Gene_id
orf_M1_M12  <- M1_M12_sig_annot$Gene_id
orf_M12_WT  <- M12_WT_sig_annot$Gene_id

# --------------------------
# Run enrichment safely
# --------------------------
res_WT      <- run_enrichment(geneid_WT, orf_WT, "WT")
res_M1      <- run_enrichment(geneid_M1, orf_M1, "M1")
res_M12     <- run_enrichment(geneid_M12, orf_M12, "M12")

res_WT_M1   <- run_enrichment(geneid_WT_M1, orf_WT_M1, "WT_M1")
res_M1_M12  <- run_enrichment(geneid_M1_M12, orf_M1_M12, "M1_M12")
res_M12_WT  <- run_enrichment(geneid_M12_WT, orf_M12_WT, "M12_WT")

# --------------------------
# Helper: save results only if not NULL
# --------------------------
save_res <- function(obj, prefix) {
  if (!is.null(obj)) {
    saveRDS(obj, file.path(OUTPUT_DIR, paste0(prefix, ".rds")))
    write.csv(as.data.frame(obj), file.path(OUTPUT_DIR, paste0(prefix, ".csv")))
  }
}

# --------------------------
# Save all results
# --------------------------
save_res(res_WT$GO,   "GO_enrichment_WT")
save_res(res_M1$GO,   "GO_enrichment_M1")
save_res(res_M12$GO,  "GO_enrichment_M12")
save_res(res_WT_M1$GO,   "GO_enrichment_WT_M1")
save_res(res_M1_M12$GO,  "GO_enrichment_M1_M12")
save_res(res_M12_WT$GO,  "GO_enrichment_M12_WT")

save_res(res_WT$KEGG,   "KEGG_enrichment_WT")
save_res(res_M1$KEGG,   "KEGG_enrichment_M1")
save_res(res_M12$KEGG,  "KEGG_enrichment_M12")
save_res(res_WT_M1$KEGG,   "KEGG_enrichment_WT_M1")
save_res(res_M1_M12$KEGG,  "KEGG_enrichment_M1_M12")
save_res(res_M12_WT$KEGG,  "KEGG_enrichment_M12_WT")

save_res(res_WT$REACT,   "Reactome_enrichment_WT")
save_res(res_M1$REACT,   "Reactome_enrichment_M1")
save_res(res_M12$REACT,  "Reactome_enrichment_M12")
save_res(res_WT_M1$REACT,   "Reactome_enrichment_WT_M1")
save_res(res_M1_M12$REACT,  "Reactome_enrichment_M1_M12")
save_res(res_M12_WT$REACT,  "Reactome_enrichment_M12_WT")

cat("Enrichment analysis completed. Results saved in:", OUTPUT_DIR, "\n")

