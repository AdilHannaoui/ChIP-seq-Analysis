# ==========================
# RIP-seq Peak Annotation (ChIPseeker)
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("config.R")

library(tidyverse)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(ggplot2)

# --------------------------
# Load TxDb from config
# --------------------------
txdb <- REFERENCE_TXDB

# --------------------------
# Load significant peak tables
# --------------------------
WT_sig_annot <- read.csv(file.path(OUTPUT_DIR, "WT_significant_peaks_with_coords.csv"))
A_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "A_significant_peaks_with_coords.csv"))
D_sig_annot  <- read.csv(file.path(OUTPUT_DIR, "D_significant_peaks_with_coords.csv"))

# --------------------------
# Convert to GRanges
# --------------------------
to_granges <- function(df) {
  GRanges(
    seqnames = df$chrom,
    ranges   = IRanges(start = df$start, end = df$end)
  )
}

gr_WT <- to_granges(WT_sig_annot)
gr_A  <- to_granges(A_sig_annot)
gr_D  <- to_granges(D_sig_annot)

# --------------------------
# Annotate peaks
# --------------------------
anno_WT <- annotatePeak(gr_WT, TxDb = txdb, tssRegion = c(-3000, 3000))
anno_A  <- annotatePeak(gr_A,  TxDb = txdb, tssRegion = c(-3000, 3000))
anno_D  <- annotatePeak(gr_D,  TxDb = txdb, tssRegion = c(-3000, 3000))

# --------------------------
# Save annotation tables
# --------------------------
write.csv(as.data.frame(anno_WT), file.path(OUTPUT_DIR, "WT_peak_annotation.csv"))
write.csv(as.data.frame(anno_A),  file.path(OUTPUT_DIR, "A_peak_annotation.csv"))
write.csv(as.data.frame(anno_D),  file.path(OUTPUT_DIR, "D_peak_annotation.csv"))

saveRDS(anno_WT, file.path(OUTPUT_DIR, "WT_peak_annotation.rds"))
saveRDS(anno_A,  file.path(OUTPUT_DIR, "A_peak_annotation.rds"))
saveRDS(anno_D,  file.path(OUTPUT_DIR, "D_peak_annotation.rds"))

# --------------------------
# Generate annotation plots
# --------------------------

# Barplot of genomic annotation categories
p_bar <- plotAnnoBar(list(WT = anno_WT, A = anno_A, D = anno_D))
ggsave(file.path(OUTPUT_DIR, "PeakAnnotation_Barplot.png"), p_bar, width = 8, height = 6)

# Pie charts for each condition
p_pie_WT <- plotAnnoPie(anno_WT)
p_pie_A  <- plotAnnoPie(anno_A)
p_pie_D  <- plotAnnoPie(anno_D)

ggsave(file.path(OUTPUT_DIR, "WT_PeakAnnotation_Pie.png"), p_pie_WT, width = 6, height = 6)
ggsave(file.path(OUTPUT_DIR, "A_PeakAnnotation_Pie.png"),  p_pie_A,  width = 6, height = 6)
ggsave(file.path(OUTPUT_DIR, "D_PeakAnnotation_Pie.png"),  p_pie_D,  width = 6, height = 6)

# Distance to TSS
p_tss <- plotDistToTSS(list(WT = anno_WT, A = anno_A, D = anno_D))
ggsave(file.path(OUTPUT_DIR, "PeakAnnotation_DistToTSS.png"), p_tss, width = 8, height = 6)

cat("Peak annotation completed. Results and plots saved in:", OUTPUT_DIR, "\n")
