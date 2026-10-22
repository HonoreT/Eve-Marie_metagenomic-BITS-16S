#!/bin/bash
#SBATCH --job-name=dada2_EMP
#SBATCH --account=def-dator3
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/dada2_EMP_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/dada2_EMP_%j.err

set -euo pipefail

module purge
module load StdEnv/2023
module load qiime2/2026.4

PROJECT=/project/166726054/Marie-Eve_Analysis
EMP=$PROJECT/04_qiime2_outputs/EMP

mkdir -p "$PROJECT/04_work_tmp"

export TMPDIR="$PROJECT/04_work_tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

echo "===== DADA2 EMP START ====="
date
hostname
echo "TMPDIR=$TMPDIR"
qiime --version

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$EMP/EMP_trimmed_demux.qza" \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 220 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-chimera-method consensus \
  --p-n-threads 16 \
  --o-table "$EMP/EMP_table.qza" \
  --o-representative-sequences "$EMP/EMP_rep-seqs.qza" \
  --o-denoising-stats "$EMP/EMP_denoising-stats.qza" \
  --o-base-transition-stats "$EMP/EMP_base-transition-stats.qza"

echo "===== DADA2 EMP END ====="
date

