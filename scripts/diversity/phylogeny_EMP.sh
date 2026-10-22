#!/bin/bash
#SBATCH --job-name=phylogeny_EMP
#SBATCH --account=def-dator3
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/phylogeny_EMP_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/phylogeny_EMP_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

PROJECT=/project/166726054/Marie-Eve_Analysis
EMP=$PROJECT/04_qiime2_outputs/EMP
PHY=$PROJECT/06_diversity/EMP/phylogeny

mkdir -p "$PHY"
mkdir -p "$PROJECT/04_work_tmp"

export TMPDIR=$PROJECT/04_work_tmp
export TMP=$TMPDIR
export TEMP=$TMPDIR

echo "===== EMP PHYLOGENY START ====="
date
hostname
qiime --version

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$EMP/EMP_rep-seqs_filtered.qza" \
  --p-n-threads 8 \
  --o-alignment "$PHY/EMP_aligned-rep-seqs.qza" \
  --o-masked-alignment "$PHY/EMP_masked-aligned-rep-seqs.qza" \
  --o-tree "$PHY/EMP_unrooted-tree.qza" \
  --o-rooted-tree "$PHY/EMP_rooted-tree.qza"

echo "===== EMP PHYLOGENY END ====="
date

