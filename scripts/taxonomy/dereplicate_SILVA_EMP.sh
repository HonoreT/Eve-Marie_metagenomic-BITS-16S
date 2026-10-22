#!/bin/bash
#SBATCH --job-name=derep_SILVA_EMP
#SBATCH --account=def-dator3
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/derep_SILVA_EMP_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/derep_SILVA_EMP_%j.err

set -euo pipefail

module load StdEnv/2023
module load qiime2/2026.4

OUT=/project/166726054/Marie-Eve_Analysis/05_taxonomy/reference_database

mkdir -p /project/166726054/Marie-Eve_Analysis/04_work_tmp
export TMPDIR=/project/166726054/Marie-Eve_Analysis/04_work_tmp
export TMP=$TMPDIR
export TEMP=$TMPDIR

echo "===== SILVA EMP DEREPLICATION START ====="
date
hostname
qiime --version

qiime rescript dereplicate \
  --i-sequences "$OUT/silva-138.2-BA-515F-806R-200-300-seqs.qza" \
  --i-taxa "$OUT/silva-138.2-BA-taxonomy.qza" \
  --p-mode lca \
  --p-perc-identity 1.0 \
  --p-threads 8 \
  --p-rank-handles domain phylum class order family genus species \
  --p-no-derep-prefix \
  --o-dereplicated-taxa "$OUT/silva-138.2-BA-515F-806R-derep-taxonomy.qza" \
  --o-dereplicated-sequences "$OUT/silva-138.2-BA-515F-806R-derep-seqs.qza"

echo "===== SILVA EMP DEREPLICATION END ====="
date

