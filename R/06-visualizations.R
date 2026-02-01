# ==========================
# RIP-seq Visualization Module
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

# --------------------------
# Libraries
# --------------------------
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(ggplot2)

# ==========================
# Load core objects
# ==========================

# DESeq2 objects
dds <- readRDS(file.path(OUTPUT_DIR, "dds.rds"))
vsd <- vst(dds, blind=FALSE)
dds_paired <- readRDS(file.path(OUTPUT_DIR, "dds_paired.rds"))
vsd_paired <- vst(dds_paired, blind=FALSE)

# DESeq2 results (contrasts)
res_list <- list(
  WT      = readRDS(file.path(OUTPUT_DIR, "DESeq2_res_WT_sig.rds")),
  M1      = readRDS(file.path(OUTPUT_DIR, "DESeq2_res_M1_sig.rds")),
  M12     = readRDS(file.path(OUTPUT_DIR, "DESeq2_res_M12_sig.rds")),
  WT_M1   = readRDS(file.path(OUTPUT_DIR, "DESeq2_res_WT_vs_M1_sig.rds")),
  M1_M12  = readRDS(file.path(OUTPUT_DIR, "DESeq2_res_M12_vs_M1_sig.rds")),
  M12_WT  = readRDS(file.path(OUTPUT_DIR, "DESeq2_res_WT_vs_M12_sig.rds"))
)

# ==========================
# Detect empty contrasts
# ==========================
valid_res <- list()

for (nm in names(res_list)) {
  res <- res_list[[nm]]
  
  if (is.null(res) || nrow(res) == 0 || !"log2FoldChange" %in% colnames(res)) {
    message("Skipping contrast ", nm, " (empty or invalid)")
  } else {
    valid_res[[nm]] <- res
  }
}

# ==========================
# Load annotation tables
# ==========================
WT_sig_annot     <- read.csv(file.path(OUTPUT_DIR, "WT_significant_peaks_with_coords.csv"))
M1_sig_annot     <- read.csv(file.path(OUTPUT_DIR, "M1_significant_peaks_with_coords.csv"))
M12_sig_annot    <- read.csv(file.path(OUTPUT_DIR, "M12_significant_peaks_with_coords.csv"))

WT_M1_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "WT_M1_significant_peaks_with_coords.csv"))
M1_M12_sig_annot <- read.csv(file.path(OUTPUT_DIR, "M1_M12_significant_peaks_with_coords.csv"))
M12_WT_sig_annot <- read.csv(file.path(OUTPUT_DIR, "M12_WT_significant_peaks_with_coords.csv"))

# ==========================
# Helper: convertir a ENTREZ
# ==========================
convert_to_entrez <- function(gene_ids) {
  if (is.null(gene_ids) || length(gene_ids) == 0) return(character(0))
  
  entrez <- mapIds(
    org.Sc.sgd.db,
    keys = as.character(gene_ids),
    column = "ENTREZID",
    keytype = GENE_ID_TYPE,
    multiVals = "first"
  )
  
  unname(entrez[!is.na(entrez)])
}

geneid_list <- list(
  WT      = convert_to_entrez(WT_sig_annot$Gene_id),
  M1      = convert_to_entrez(M1_sig_annot$Gene_id),
  M12     = convert_to_entrez(M12_sig_annot$Gene_id),
  WT_M1   = convert_to_entrez(WT_M1_sig_annot$Gene_id),
  M1_M12  = convert_to_entrez(M1_M12_sig_annot$Gene_id),
  M12_WT  = convert_to_entrez(M12_WT_sig_annot$Gene_id)
)

# ==========================
# PCA
# ==========================
ggsave(file.path(OUTPUT_DIR, "PCA_plot.png"),
       plotPCA(vsd, intgroup="condition") + ggtitle("PCA of RIP-seq samples"),
       width=7, height=6)

ggsave(file.path(OUTPUT_DIR, "PCA_paired_plot.png"),
       plotPCA(vsd_paired, intgroup="condition") + ggtitle("PCA of RIP-seq samples"),
       width=7, height=6)

# ==========================
# MA plots (only valid contrasts)
# ==========================
for (nm in names(valid_res)) {
  png(file.path(OUTPUT_DIR, paste0("MAplot_", nm, ".png")), width=1200, height=900)
  plotMA(valid_res[[nm]], ylim=c(-5,5), cex=1.2)
  dev.off()
}

# ==========================
# Correlation heatmaps
# ==========================
png(file.path(OUTPUT_DIR, "CorrelationHeatmap_samples.png"), width=1200, height=1000)
Heatmap(cor(assay(vsd)), name="correlation")
dev.off()

png(file.path(OUTPUT_DIR, "CorrelationHeatmap_paired_samples.png"), width=1200, height=1000)
Heatmap(cor(assay(vsd_paired)), name="correlation")
dev.off()

# ==========================
# UpSet plots (only valid contrasts)
# ==========================
sig_genes <- function(res) {
  if (is.null(res)) return(character(0))
  rownames(res[which(res$padj < 0.05 & !is.na(res$padj)), ])
}

# Unpaired
unpaired_valid <- intersect(c("WT","M1","M12"), names(valid_res))

if (length(unpaired_valid) > 0) {
  upset_list <- lapply(unpaired_valid, function(nm) sig_genes(valid_res[[nm]]))
  names(upset_list) <- unpaired_valid
  
  png(file.path(OUTPUT_DIR, "UpSet_WT_M1_M12.png"), width=1200, height=900)
  upset(fromList(upset_list), order.by="freq")
  dev.off()
}

# Paired
paired_valid <- intersect(c("WT_M1","M1_M12","M12_WT"), names(valid_res))

if (length(paired_valid) > 0) {
  upset_list_paired <- lapply(paired_valid, function(nm) sig_genes(valid_res[[nm]]))
  names(upset_list_paired) <- paired_valid
  
  png(file.path(OUTPUT_DIR, "UpSet_paired_WT_M1_M12.png"), width=1200, height=900)
  upset(fromList(upset_list_paired), order.by="freq")
  dev.off()
}

# ==========================
# SAFE DOTPLOT FUNCTION
# ==========================
safe_dotplot <- function(geneid_subset, filename, fun, ...) {
  # Filtrar grupos vacíos
  geneid_subset <- geneid_subset[sapply(geneid_subset, function(x) length(x) > 0)]
  
  if (length(geneid_subset) == 0) {
    message("Skipping ", filename, " (all gene lists empty)")
    return()
  }
  
  # Intentar enriquecimiento
  cc <- try(compareCluster(geneid_subset, fun = fun, ...), silent = TRUE)
  
  if (inherits(cc, "try-error") || is.null(cc) || nrow(cc) == 0) {
    message("Skipping ", filename, " (no enrichment found)")
    return()
  }
  
  # Intentar dotplot
  dp <- try(dotplot(cc), silent = TRUE)
  
  if (inherits(dp, "try-error")) {
    message("Skipping ", filename, " (dotplot failed)")
    return()
  }
  
  ggsave(file.path(OUTPUT_DIR, filename), dp, width = 10, height = 8)
}

# ==========================
# Dotplots (GO, KEGG, Reactome)
# ==========================

# GO
safe_dotplot(
  geneid_list[c("WT","M1","M12")],
  "Dotplot_GO.png",
  fun = "enrichGO",
  OrgDb = org.Sc.sgd.db,
  ont = GO_ONTOLOGY
)

# GO paired
safe_dotplot(
  geneid_list[c("WT_M1","M1_M12","M12_WT")],
  "Dotplot_GO_paired.png",
  fun = "enrichGO",
  OrgDb = org.Sc.sgd.db,
  ont = GO_ONTOLOGY
)

# KEGG
safe_dotplot(
  geneid_list[c("WT","M1","M12")],
  "Dotplot_KEGG.png",
  fun = "enrichKEGG",
  organism = KEGG_ORGANISM
)

# KEGG paired
safe_dotplot(
  geneid_list[c("WT_M1","M1_M12","M12_WT")],
  "Dotplot_KEGG_paired.png",
  fun = "enrichKEGG",
  organism = KEGG_ORGANISM
)

# Reactome
safe_dotplot(
  geneid_list[c("WT","M1","M12")],
  "Dotplot_Reactome.png",
  fun = "enrichPathway",
  organism = REACTOME_ORGANISM
)

# Reactome paired
safe_dotplot(
  geneid_list[c("WT_M1","M1_M12","M12_WT")],
  "Dotplot_Reactome_paired.png",
  fun = "enrichPathway",
  organism = REACTOME_ORGANISM
)

cat("Visualization module completed. All plots saved in:", OUTPUT_DIR, "\n")
