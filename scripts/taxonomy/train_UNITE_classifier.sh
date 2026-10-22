#!/bin/bash
#SBATCH --job-name=UNITE_NB
#SBATCH --account=def-dator3
##SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/UNITE_NB_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/UNITE_NB_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

PROJECT=/project/166726054/Marie-Eve_Analysis
OUT="$PROJECT/05_taxonomy/UNITE_2025-02-19_dynamic"

TMPBASE="$PROJECT/tmp_qiime"
mkdir -p "$TMPBASE"
export TMPDIR="$TMPBASE"

CACHE="$PROJECT/tmp_qiime/qiime_cache_UNITE"
mkdir -p "$CACHE"

REFSEQ="$OUT/unite-2025-02-19-eukaryotes-dynamic-seqs-QIIME2-2026.4.qza"
REFTAX="$OUT/unite-2025-02-19-eukaryotes-dynamic-taxonomy-QIIME2-2026.4.qza"

CLASSIFIER="$OUT/unite-2025-02-19-eukaryotes-dynamic-classifier-QIIME2-2026.4.qza"

echo "===== Environment ====="
echo "HOSTNAME=$(hostname)"
echo "TMPDIR=$TMPDIR"
echo "REFSEQ=$REFSEQ"
echo "REFTAX=$REFTAX"
echo

echo "===== Input artifacts ====="
qiime tools peek "$REFSEQ"
qiime tools peek "$REFTAX"

echo
echo "===== Training classifier ====="

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "$REFSEQ" \
  --i-reference-taxonomy "$REFTAX" \
  --o-classifier "$CLASSIFIER" \
  --use-cache "$CACHE"

echo
echo "===== Classifier created ====="
qiime tools peek "$CLASSIFIER"
ls -lh "$CLASSIFIER"


