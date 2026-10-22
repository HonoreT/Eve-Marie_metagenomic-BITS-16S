#!/bin/bash
#SBATCH --job-name=classify_EMP
#SBATCH --account=def-dator3
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/classify_EMP_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/classify_EMP_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

PROJECT=/project/166726054/Marie-Eve_Analysis
EMP=$PROJECT/04_qiime2_outputs/EMP
REF=$PROJECT/05_taxonomy/reference_database
TAX=$PROJECT/05_taxonomy

mkdir -p "$PROJECT/04_work_tmp"
export TMPDIR="$PROJECT/04_work_tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

echo "===== EMP TAXONOMY CLASSIFICATION START ====="
date
hostname
qiime --version

qiime feature-classifier classify-sklearn \
  --i-classifier "$REF/silva-138.2-BA-515F-806R-classifier-QIIME2-2026.4.qza" \
  --i-reads "$EMP/EMP_rep-seqs.qza" \
  --p-n-jobs 8 \
  --o-classification "$TAX/EMP_taxonomy.qza"

echo "===== EMP TAXONOMY CLASSIFICATION END ====="
date

