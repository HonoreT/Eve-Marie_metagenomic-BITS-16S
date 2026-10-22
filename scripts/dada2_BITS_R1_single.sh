#!/bin/bash
#SBATCH --job-name=BITS_R1_DADA2
#SBATCH --account=def-dator3
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/BITS_R1_DADA2_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/BITS_R1_DADA2_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

PROJECT=/project/166726054/Marie-Eve_Analysis
Q2BITS="$PROJECT/04_qiime2_outputs/BITS"

qiime dada2 denoise-single \
  --i-demultiplexed-seqs "$Q2BITS/BITS_R1_trimmed.qza" \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --p-max-ee 2 \
  --p-trunc-q 2 \
  --p-chimera-method consensus \
  --p-n-threads 8 \
  --o-table "$Q2BITS/BITS_table_R1_single_maxEE2.qza" \
  --o-representative-sequences "$Q2BITS/BITS_rep_seqs_R1_single_maxEE2.qza" \
  --o-denoising-stats "$Q2BITS/BITS_denoising_stats_R1_single_maxEE2.qza" \
  --o-base-transition-stats "$Q2BITS/BITS_base_transition_stats_R1_single_maxEE2.qza"

