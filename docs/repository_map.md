# Repository map

## Notebooks

`notebooks/`  
Google Colab / Jupyter notebooks used for exploratory and reproducible analysis.

The original `Eve-Marie_Analyses.ipynb` notebook should be saved here from Google Colab. It is intentionally not reconstructed from screenshots or partial records.

## Scripts

### 16S / EMP
`scripts/16S_EMP/`

- DADA2 denoising for bacterial/archaeal amplicons.

### ITS / BITS
`scripts/ITS_BITS/`

- post-primer quality control
- paired-end diagnostic workflows
- final R1 single-end DADA2 workflow

### Taxonomy
`scripts/taxonomy/`

- SILVA 138.2 preparation/training/classification
- UNITE 2025-02-19 training/classification

### Diversity
`scripts/diversity/`

- phylogeny
- alpha rarefaction
- core diversity metrics
- ITS diversity analyses

### Differential abundance
`scripts/differential_abundance/`

- ANCOM-BC2 at ASV and genus levels

### PICRUSt2 / MetaCyc
`scripts/PICRUSt2_MetaCyc/`

- PICRUSt2 workflow
- MetaCyc functional summaries
- targeted limma analyses

### FUNGuild
`scripts/FUNGuild/`

- ecological guild inference
- trophic-mode figures
- targeted limma analyses
- HTML outputs and volcano plots

## Metadata

`metadata/`

Only non-sensitive metadata, templates, and data dictionaries should be committed.

## Figures

`figures/`

Final lightweight publication or presentation figures may be added here. Large interactive HTML files should normally remain outside the repository unless specifically needed.

## Data not stored in GitHub

The following remain outside GitHub:

- raw FASTQ files
- large QIIME 2 artifacts (`.qza`, `.qzv`)
- temporary working directories
- large intermediate PICRUSt2 outputs
- confidential or identifying metadata
