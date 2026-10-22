#!/bin/bash
#SBATCH --job-name=FastQC_BITS
#SBATCH --account=def-dator3
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/FastQC_BITS_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/FastQC_BITS_%j.err

set -euo pipefail

module load StdEnv/2023
module load fastqc/0.12.1

PROJECT=/project/166726054/Marie-Eve_Analysis
EXPORT=$PROJECT/04_qiime2_outputs/BITS/BITS_trimmed_export
FASTQC=$PROJECT/03_quality_control/BITS_post_primer/fastqc

mkdir -p "$FASTQC"

echo "Nombre de FASTQ :"
find "$EXPORT" -type f -name "*.fastq.gz" | wc -l

fastqc \
  --threads 8 \
  --outdir "$FASTQC" \
  "$EXPORT"/*.fastq.gz

echo "FastQC terminé"
date

find "$FASTQC" -type f | wc -l

