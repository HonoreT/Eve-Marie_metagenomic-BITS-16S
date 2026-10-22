#!/bin/bash
#SBATCH --job-name=PICRUSt2_EMP
#SBATCH --account=def-dator3
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/project/166726054/Marie-Eve_Analysis/12_logs/PICRUSt2_EMP_%j.out
#SBATCH --error=/project/166726054/Marie-Eve_Analysis/12_logs/PICRUSt2_EMP_%j.err

set -euo pipefail

module load StdEnv/2023
module load apptainer/1.3.5

PROJECT=/project/166726054/Marie-Eve_Analysis
FUNC=$PROJECT/08_functional_analysis/EMP

SIF=$FUNC/container/picrust2_2.6.3-2.sif
FASTA=$FUNC/input/rep-seqs-export/dna-sequences.fasta
BIOM=$FUNC/input/table-export/feature-table.biom
OUT=$FUNC/picrust2_output

mkdir -p "$FUNC/tmp"
mkdir -p "$FUNC/container/cache"

export APPTAINER_CACHEDIR="$FUNC/container/cache"
export APPTAINER_TMPDIR="$FUNC/tmp"

# PICRUSt2 préfère créer son répertoire de sortie lui-même.
rm -rf "$OUT"

echo "===== PICRUSt2 EMP START ====="
date
hostname
apptainer --version

echo "FASTA: $FASTA"
echo "BIOM : $BIOM"
echo "OUT  : $OUT"

apptainer exec \
  --bind "$PROJECT:$PROJECT" \
  "$SIF" \
  picrust2_pipeline.py \
    -s "$FASTA" \
    -i "$BIOM" \
    -o "$OUT" \
    -p 8 \
    --verbose

echo "===== PICRUSt2 EMP END ====="
date

