# Eve-Marie_metagenomic-BITS-16S

Reproducible bioinformatics workflows for the Marie-Ève soil microbiome project.

## Scope

This repository documents the analysis of two amplicon datasets:

- **EMP / 16S**: bacterial and archaeal communities
- **BITS / ITS**: fungal communities

It also documents downstream functional/ecological inference using:

- **PICRUSt2 / MetaCyc** for predicted prokaryotic functional potential
- **FUNGuild** for taxonomically inferred fungal ecological guilds

## Experimental design

- 48 soil samples
- 3 cultivars: Blueray, Duke, Patriot
- 2 fertilization regimes: Mineral, Organic
- 2 PGPR conditions: No PGPR, PGPR
- 4 treatments: CTRLMin, PGPRMin, CTROrg, PGPROrg
- 4 blocks

Statistical framework: `Fertilisation × PGPR × Cultivar + Bloc`.

## Repository structure

```text
.
├── README.md
├── .gitignore
├── notebooks/
├── scripts/
│   ├── 16S_EMP/
│   ├── ITS_BITS/
│   ├── taxonomy/
│   ├── diversity/
│   ├── differential_abundance/
│   ├── PICRUSt2_MetaCyc/
│   └── FUNGuild/
├── metadata/
├── figures/
├── docs/
└── environment/
```

## Important data policy

Raw sequencing data and large intermediate QIIME 2 files are intentionally **not stored in this repository**. In particular, do not commit:

- raw FASTQ files
- `.qza` / `.qzv` files unless specifically needed and size-appropriate
- temporary working directories
- large BIOM tables or intermediate PICRUSt2 outputs
- credentials, tokens, or sensitive metadata

## Current analysis environment

- QIIME 2: 2026.4
- PICRUSt2: 2.6.3
- Apptainer: 1.3.5
- Reference databases: SILVA 138.2 and UNITE 2025-02-19

See `docs/workflow_methodology.md` and `environment/software_versions.txt` for details.

## Status

The repository structure and documentation have been initialized. The validated scripts from the Valeria project directory `11_scripts/` should be copied into the corresponding subdirectories next.
