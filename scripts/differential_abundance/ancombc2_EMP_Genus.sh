#!/bin/bash
#SBATCH --job-name=abc2_EMP_GEN
#SBATCH --account=def-dator3
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/abc2_EMP_GEN_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/abc2_EMP_GEN_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

PROJECT=/project/166726054/Marie-Eve_Analysis
DA=$PROJECT/07_differential_abundance/EMP
GENUS=$DA/Genus

TABLE=$GENUS/EMP_genus_table.qza
META=$DA/Metadata_ANCOMBC2_prefixed.tsv

mkdir -p "$PROJECT/04_work_tmp"

export TMPDIR=$PROJECT/04_work_tmp
export TMP=$TMPDIR
export TEMP=$TMPDIR

echo "===== ANCOM-BC2 EMP GENUS START ====="
date
hostname
qiime --version

qiime composition ancombc2 \
  --i-table "$TABLE" \
  --m-metadata-file "$META" \
  --p-fixed-effects-formula 'Traitement * Cultivar + Bloc' \
  --p-reference-levels \
      Traitement::CTRLMin \
      Cultivar::Blueray \
      Bloc::1 \
  --p-p-adjust-method BH \
  --p-prevalence-cutoff 0.10 \
  --p-alpha 0.05 \
  --p-diff-robust \
  --p-num-processes 8 \
  --o-ancombc2-output "$GENUS/EMP_Genus_ANCOMBC2.qza"

echo "===== ANCOM-BC2 EMP GENUS END ====="
date

