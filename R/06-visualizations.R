# ==========================
# RIP-seq Visualization Module
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("config.R")

# --------------------------
# Libraries
# --------------------------
library(DESeq2)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(tidyverse)
library(ggplot2)

# ==========================
# Load core objects
# ==========================

# DESeq2 objects
dds <- readRDS(file.path(OUTPUT_DIR, "dds_object.rds"))
vsd <- readRDS(file.path(OUTPUT_DIR, "vsd_object.rds"))

# DESeq2 results (contrasts)
res_WT <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_WT.rds"))  # si lo usas
res_A  <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_A.rds"))
res_D  <- readRDS(file.path(OUTPUT_DIR, "DESeq2_res_D.rds"))

# Tablas de picos significativos (para listas de genes)
WT_sig_annot <- read.csv(file.path(OUTPUT_DIR, "WT_significant_peaks_with_coords.csv"))
A_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "A_significant_peaks_with_coords.csv"))
D_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "D_significant_peaks_with_coords.csv"))

# ==========================
# Helper: convertir a ENTREZ
# ==========================
convert_to_entrez <- function(gene_ids) {
  entrez <- mapIds(
    org.Sc.sgd.db,
    keys     = as.character(gene_ids),
    column   = ENTREZ_COLUMN,   # definido en config.R, normalmente "ENTREZID"
    keytype  = GENE_ID_TYPE,    # definido en config.R, p.ej. "ENSEMBL" o "SYMBOL"
    multiVals = "first"
  )
  unname(entrez[!is.na(entrez)])
}

geneid_WT <- convert_to_entrez(WT_sig_annot$Gene_id)
geneid_A  <- convert_to_entrez(A_sig_annot$Gene_id)
geneid_D  <- convert_to_entrez(D_sig_annot$Gene_id)

# ==========================
# PCA
# ==========================
pca <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA of RIP-seq samples")

ggsave(file.path(OUTPUT_DIR, "PCA_plot.png"), pca, width = 7, height = 6)

# ==========================
# MA plots (WT vs A, WT vs D)
# ==========================
png(file.path(OUTPUT_DIR, "MAplot_WT_vs_A.png"), width = 1200, height = 900)
plotMA(res_A, ylim = c(-5, 5), main = "MA Plot WT vs A")
dev.off()

png(file.path(OUTPUT_DIR, "MAplot_WT_vs_D.png"), width = 1200, height = 900)
plotMA(res_D, ylim = c(-5, 5), main = "MA Plot WT vs D")
dev.off()

# ==========================
# Enhanced Volcano plots
# ==========================
volcano_A <- EnhancedVolcano(
  res_A,
  lab = rownames(res_A),
  x = "log2FoldChange",
  y = "padj",
  title = "Volcano Plot WT vs A"
)

volcano_D <- EnhancedVolcano(
  res_D,
  lab = rownames(res_D),
  x = "log2FoldChange",
  y = "padj",
  title = "Volcano Plot WT vs D"
)

ggsave(file.path(OUTPUT_DIR, "Volcano_WT_vs_A.png"), volcano_A, width = 8, height = 7)
ggsave(file.path(OUTPUT_DIR, "Volcano_WT_vs_D.png"), volcano_D, width = 8, height = 7)

# ==========================
# Sample-to-sample correlation heatmap
# ==========================
vsd_mat <- assay(vsd)
cor_mat <- cor(vsd_mat)

png(file.path(OUTPUT_DIR, "CorrelationHeatmap_samples.png"), width = 1200, height = 1000)
ht <- Heatmap(
  cor_mat,
  name = "correlation",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Sample-to-sample correlation",
  row_title = "Samples",
  col = colorRamp2(c(0.7, 0.85, 1), c("white", "skyblue", "darkblue"))
)
draw(ht)
dev.off()

# ==========================
# UpSet plot (intersección de genes significativos)
# ==========================
WT_sig_genes <- rownames(res_WT[which(res_WT$padj < 0.05 & !is.na(res_WT$padj)), ])
A_sig_genes  <- rownames(res_A[which(res_A$padj < 0.05 & !is.na(res_A$padj)), ])
D_sig_genes  <- rownames(res_D[which(res_D$padj < 0.05 & !is.na(res_D$padj)), ])

upset_list <- list(
  WT = WT_sig_genes,
  A  = A_sig_genes,
  D  = D_sig_genes
)

png(file.path(OUTPUT_DIR, "UpSet_WT_A_D.png"), width = 1200, height = 900)
upset(fromList(upset_list), order.by = "freq", main.bar.color = "steelblue")
dev.off()

# ==========================
# Dotplots de enriquecimiento (GO, KEGG, Reactome)
# ==========================
# GO
dot_GO <- dotplot(
  compareCluster(
    list(WT = geneid_WT, A = geneid_A, D = geneid_D),
    fun   = "enrichGO",
    OrgDb = org.Sc.sgd.db,
    ont   = GO_ONTOLOGY
  )
) + ggtitle("GO enrichment comparison (WT, A, D)")

ggsave(file.path(OUTPUT_DIR, "Dotplot_GO.png"), dot_GO, width = 10, height = 8)

# KEGG
dot_KEGG <- dotplot(
  compareCluster(
    list(WT = geneid_WT, A = geneid_A, D = geneid_D),
    fun      = "enrichKEGG",
    organism = KEGG_ORGANISM
  )
) + ggtitle("KEGG enrichment comparison (WT, A, D)")

ggsave(file.path(OUTPUT_DIR, "Dotplot_KEGG.png"), dot_KEGG, width = 10, height = 8)

# Reactome
dot_Reactome <- dotplot(
  compareCluster(
    list(WT = geneid_WT, A = geneid_A, D = geneid_D),
    fun      = "enrichPathway",
    organism = REACTOME_ORGANISM
  )
) + ggtitle("Reactome enrichment comparison (WT, A, D)")

ggsave(file.path(OUTPUT_DIR, "Dotplot_Reactome.png"), dot_Reactome, width = 10, height = 8)

cat("Visualization module completed. All plots saved in:", OUTPUT_DIR, "\n")
