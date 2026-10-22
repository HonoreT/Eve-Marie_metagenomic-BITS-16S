#!/bin/bash
#SBATCH --job-name=DADA2_BITS
#SBATCH --account=def-dator3
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/DADA2_BITS_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/DADA2_BITS_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

PROJECT=/project/166726054/Marie-Eve_Analysis
Q2BITS=$PROJECT/04_qiime2_outputs/BITS
TMPQ=$PROJECT/04_work_tmp

mkdir -p "$TMPQ"

export TMPDIR="$TMPQ"
export TMP="$TMPQ"
export TEMP="$TMPQ"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$Q2BITS/BITS_primers_trimmed.qza" \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-chimera-method consensus \
  --p-n-threads 16 \
  --o-table "$Q2BITS/BITS_table_full_length.qza" \
  --o-representative-sequences "$Q2BITS/BITS_rep_seqs_full_length.qza" \
  --o-denoising-stats "$Q2BITS/BITS_denoising_stats_full_length.qza" \
  --o-base-transition-stats "$Q2BITS/BITS_base_transition_stats_full_length.qza"

qiime metadata tabulate \
  --m-input-file "$Q2BITS/BITS_denoising_stats_full_length.qza" \
  --o-visualization "$Q2BITS/BITS_denoising_stats_full_length.qzv"

