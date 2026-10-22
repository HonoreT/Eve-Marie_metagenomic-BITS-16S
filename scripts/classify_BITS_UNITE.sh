#!/bin/bash
#SBATCH --job-name=BITS_UNITE
#SBATCH --account=def-dator3
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/BITS_UNITE_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/BITS_UNITE_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

PROJECT=/project/166726054/Marie-Eve_Analysis

Q2BITS="$PROJECT/04_qiime2_outputs/BITS"
TAXOUT="$PROJECT/05_taxonomy/BITS_UNITE"

CLASSIFIER="$PROJECT/05_taxonomy/UNITE_2025-02-19_dynamic/unite-2025-02-19-eukaryotes-dynamic-classifier-QIIME2-2026.4.qza"
REPSEQS="$Q2BITS/BITS_rep_seqs_R1_single_maxEE2.qza"

TMPBASE="$PROJECT/tmp_qiime"
mkdir -p "$TMPBASE"
export TMPDIR="$TMPBASE"

CACHE="$PROJECT/tmp_qiime/qiime_cache_BITS_UNITE"

mkdir -p "$TAXOUT"

echo "===== INPUTS ====="
echo "HOSTNAME=$(hostname)"
echo "CLASSIFIER=$CLASSIFIER"
echo "REPSEQS=$REPSEQS"
echo "TAXOUT=$TAXOUT"
echo "TMPDIR=$TMPDIR"
echo "CACHE=$CACHE"
echo

echo "===== CHECK ARTIFACTS ====="
qiime tools peek "$CLASSIFIER"
qiime tools peek "$REPSEQS"

echo
echo "===== CLASSIFICATION ====="

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$REPSEQS" \
  --p-n-jobs 1 \
  --p-reads-per-batch 800 \
  --p-confidence 0.7 \
  --p-read-orientation auto \
  --o-classification "$TAXOUT/BITS_taxonomy_UNITE.qza" \
  --use-cache "$CACHE"

echo
echo "===== TAXONOMY VISUALIZATION ====="

qiime metadata tabulate \
  --m-input-file "$TAXOUT/BITS_taxonomy_UNITE.qza" \
  --o-visualization "$TAXOUT/BITS_taxonomy_UNITE.qzv"

echo
echo "===== EXPORT TAXONOMY ====="

rm -rf "$TAXOUT/exported"

qiime tools export \
  --input-path "$TAXOUT/BITS_taxonomy_UNITE.qza" \
  --output-path "$TAXOUT/exported"

echo
echo "===== COMPLETE ====="

qiime tools peek "$TAXOUT/BITS_taxonomy_UNITE.qza"

ls -lh \
  "$TAXOUT/BITS_taxonomy_UNITE.qza" \
  "$TAXOUT/BITS_taxonomy_UNITE.qzv" \
  "$TAXOUT/exported/taxonomy.tsv"
