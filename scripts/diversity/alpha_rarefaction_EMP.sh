#!/bin/bash
#SBATCH --job-name=alpha_rare_EMP
#SBATCH --account=def-dator3
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/alpha_rare_EMP_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/alpha_rare_EMP_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

PROJECT=/project/166726054/Marie-Eve_Analysis
EMP=$PROJECT/04_qiime2_outputs/EMP
DIV=$PROJECT/06_diversity/EMP

mkdir -p "$DIV"
mkdir -p "$PROJECT/04_work_tmp"

export TMPDIR="$PROJECT/04_work_tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

echo "===== EMP ALPHA RAREFACTION START ====="
date
hostname
qiime --version

qiime diversity alpha-rarefaction \
  --i-table "$EMP/EMP_table_filtered.qza" \
  --p-min-depth 1000 \
  --p-max-depth 30000 \
  --p-metrics observed_features shannon \
  --p-steps 20 \
  --p-iterations 10 \
  --o-visualization "$DIV/EMP_alpha-rarefaction_30000.qzv"

echo "===== EMP ALPHA RAREFACTION END ====="
date

