#!/bin/bash

# ==================================================
# Global configuration file for RIP-seq pipeline
# ==================================================

# ----------------------
# Project directories
# ----------------------
WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTQ_DIR="$WORKDIR/data"
OUTPUT_DIR="$WORKDIR/output"

# ----------------------
# Computational resources
# ----------------------
THREADS=6

# ----------------------
# Tools and references
# ----------------------

# Trimmomatic
TRIMMO_JAR="$WORKDIR/Trimmomatic-0.39/trimmomatic-0.39.jar"
TRIMMOMATIC_ADAPTERS="$WORKDIR/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"

# BOWTIE2
HISAT2_INDEX="$WORKDIR/BOWTIE2/cerevisiae/index/genome"

# ----------------------
# Output subdirectories
# ----------------------
FASTQC_PRE_DIR="$OUTPUT_DIR/fastqc_pre"
FASTQC_POST_DIR="$OUTPUT_DIR/fastqc_post"
TRIMMED_FASTQ_DIR="$OUTPUT_DIR/fastq_trimmed"
BOWTIE2_DIR="$OUTPUT_DIR/bowtie2"
LOG_DIR="$OUTPUT_DIR/logs"
