#!/bin/bash
#SBATCH --job-name=core_EMP
#SBATCH --account=def-dator3
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/core_EMP_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/core_EMP_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

PROJECT=/project/166726054/Marie-Eve_Analysis
FINAL=$PROJECT/06_exports/EMP_final
OUT=$PROJECT/06_diversity/EMP/core_metrics_30000

mkdir -p "$PROJECT/04_work_tmp"

# QIIME 2 exige que --output-dir n'existe pas déjà
rm -rf "$OUT"

export TMPDIR=$PROJECT/04_work_tmp
export TMP=$TMPDIR
export TEMP=$TMPDIR

echo "===== EMP CORE METRICS START ====="
date
hostname
qiime --version

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$FINAL/rooted-tree.qza" \
  --i-table "$FINAL/table-emp-prok.qza" \
  --p-sampling-depth 30000 \
  --m-metadata-file "$FINAL/Metadata.tsv" \
  --p-n-jobs-or-threads 8 \
  --output-dir "$OUT"

echo "===== EMP CORE METRICS END ====="
date

