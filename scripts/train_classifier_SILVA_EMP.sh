#!/bin/bash
#SBATCH --job-name=train_SILVA_EMP
#SBATCH --account=def-dator3
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/train_SILVA_EMP_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/train_SILVA_EMP_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

OUT=/project/166726054/Marie-Eve_Analysis/05_taxonomy/reference_database

mkdir -p /project/166726054/Marie-Eve_Analysis/04_work_tmp
export TMPDIR=/project/166726054/Marie-Eve_Analysis/04_work_tmp
export TMP=$TMPDIR
export TEMP=$TMPDIR

echo "===== SILVA EMP CLASSIFIER TRAINING START ====="
date
hostname
qiime --version

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "$OUT/silva-138.2-BA-515F-806R-derep-seqs.qza" \
  --i-reference-taxonomy "$OUT/silva-138.2-BA-515F-806R-derep-taxonomy.qza" \
  --o-classifier "$OUT/silva-138.2-BA-515F-806R-classifier-QIIME2-2026.4.qza"

echo "===== SILVA EMP CLASSIFIER TRAINING END ====="
date


