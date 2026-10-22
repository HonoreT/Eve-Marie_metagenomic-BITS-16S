#!/bin/bash
#SBATCH --job-name=BITS_DIV
#SBATCH --account=def-dator3
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/BITS_DIV_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/BITS_DIV_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

PROJECT=/project/166726054/Marie-Eve_Analysis
BITSFINAL="$PROJECT/06_exports/BITS_final"
DIVBITS="$PROJECT/06_diversity/BITS"

TMPBASE="$PROJECT/tmp_qiime"
mkdir -p "$TMPBASE"
export TMPDIR="$TMPBASE"

mkdir -p "$DIVBITS"

OUT="$DIVBITS/core_metrics_12000"

# QIIME2 --output-dir exige que le dossier cible n'existe pas déjà
rm -rf "$OUT"

echo "===== INPUT ====="
qiime tools peek "$BITSFINAL/table-bits-fungi.qza"

echo
echo "===== CORE METRICS: sampling depth 12000 ====="

qiime diversity core-metrics \
  --i-table "$BITSFINAL/table-bits-fungi.qza" \
  --p-sampling-depth 12000 \
  --m-metadata-file "$BITSFINAL/Metadata.tsv" \
  --p-n-jobs 4 \
  --output-dir "$OUT"

echo
echo "===== COMPLETE ====="

find "$OUT" -maxdepth 1 -type f -printf '%f\n' | sort
