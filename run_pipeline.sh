#!/usr/bin/env bash

set -e
set -o pipefail

SECONDS=0
echo "=== Pipeline started successfully ===)"

# 1. Environment Activation
echo "[1/10] Activating Conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ripseq-rpb4

# 2. Configuration loading
echo "[2/10] Loading R config..."
CONFIG_R="R/config.R"
echo "[2/7] Loading bash config..."
CONFIG_bash="bash/config.sh"


# 3. QC
echo "[3/10] Running FastQC..."
bash bash/01-fastqc.sh $CONFIG_bash

# 4. Trimming
echo "[4/10] Running Trimmomatic..."
bash bash/02-trimming.sh $CONFIG_bash

# 5. Alignment
echo "[5/10] Running Bowtie2..."
bash bash/03-alignment_bowtie2.sh $CONFIG_bash

# 6. Peak Calling
echo "[6/10] Running macs2..."
bash bash/04-peak_calling.sh $CONFIG_bash

# 7. Peak Intersection
echo "[7/10] Making Peak Intersection..."
bash bash/05-Peak_intersection.sh $CONFIG_bash

# 8. R Analysis
echo "[8/10] Running DESeq2 and downstream analysis..."
Rscript R/01-deseq2_analysis.R $CONFIG_R
Rscript R/02-coordinates_extraction.R $CONFIG_R
Rscript R/03-Sequence_extraction.R $CONFIG_R
Rscript R/04-enrichment_analysis.R $CONFIG_R
Rscript R/05-peak_annotation.R $CONFIG_R
Rscript R/06-visualizations.R $CONFIG_R

# 9. Motif Analysis
echo "[9/10] Finding Motifs..."
bash bash/06-Motif_analysis.sh $CONFIG_bash

# 10. IGV Running
echo "[10/10] Running IGV..."
bash bash/07-igv_running.sh $CONFIG_bash


duration=$SECONDS
echo "Total execution time: $((duration / 60)) minutes $((duration % 60)) seconds"
