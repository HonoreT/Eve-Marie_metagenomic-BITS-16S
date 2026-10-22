# Workflow methodology

## Project overview

The project analyzes soil microbiome data from 48 samples using two marker-gene workflows:

1. **EMP / 16S** for bacterial and archaeal communities.
2. **BITS / ITS** for fungal communities.

The downstream analyses include taxonomic profiling, alpha and beta diversity, differential abundance, predicted prokaryotic functional potential with PICRUSt2/MetaCyc, and fungal ecological guild inference with FUNGuild.

## Experimental design

- Cultivars: Blueray, Duke, Patriot
- Fertilization: Mineral, Organic
- PGPR: No PGPR, PGPR
- Treatments: CTRLMin, PGPRMin, CTROrg, PGPROrg
- Blocks: 4

Primary statistical framework:

```text
Fertilisation × PGPR × Cultivar + Bloc
```

## EMP / 16S

### Denoising

- QIIME 2: 2026.4
- DADA2 paired-end
- primers: 515F / 806R
- trim-left-f: 0
- trim-left-r: 0
- trunc-len-f: 250
- trunc-len-r: 220
- max-ee-f: 2
- max-ee-r: 2
- chimera method: consensus
- final dataset: 48 samples, 7,553 ASVs

### Taxonomy

- SILVA 138.2
- region-specific classifier for 515F/806R
- Bacteria + Archaea

### Diversity

Rarefaction depth: 30,000 reads/sample.

Alpha diversity:
- Observed features
- Shannon
- Pielou evenness
- Faith's PD

Beta diversity:
- Bray-Curtis
- Jaccard
- unweighted UniFrac
- weighted UniFrac

### Differential abundance

ANCOM-BC2 was used at ASV and genus levels.

Working model:

```text
Treatment × Cultivar + Bloc
```

Reference levels:
- CTRLMin
- Blueray
- Bloc 1

Multiple-testing correction:
- Benjamini-Hochberg FDR

## PICRUSt2 / MetaCyc

- PICRUSt2: 2.6.3
- Apptainer: 1.3.5
- CPUs: 8

Conceptual command:

```bash
picrust2_pipeline.py -s <FASTA> -i <BIOM> -o <OUT> -p 8 --verbose
```

Outputs:
- 2,995 EC
- 9,105 KO
- 563 MetaCyc pathways

After prevalence filtering at >=25%:
- 519 MetaCyc pathways retained

Targeted functional signatures:
- denitrification-related
- assimilatory nitrate reduction
- allantoin degradation
- phosphonate degradation
- aerobactin biosynthesis
- 2,3-butanediol biosynthesis

Targeted statistical model:

```r
~ Fertilisation * PGPR + Cultivar + Bloc
```

Additive model for main effects:

```r
~ Fertilisation + PGPR + Cultivar + Bloc
```

## BITS / ITS

Paired-end processing was evaluated first, but the reverse-read quality/retention was inadequate. The final workflow therefore used R1 single-end reads.

Final DADA2 settings:
- single-end R1
- trim-left: 0
- trunc-len: 0
- max-ee: 2
- chimera method: consensus

Final output:
- 48 samples
- 1,034,406 non-chimeric reads
- 6,381 ASVs

### Taxonomy

- UNITE 2025-02-19
- dynamic taxonomy
- eukaryotes

Fungi-only dataset:
- 5,272 fungal ASVs
- 986,462 fungal reads
- 48 samples
- minimum sample depth: 12,261 reads

### ITS diversity

Rarefaction depth: 12,000 reads/sample.

Alpha diversity:
- Observed features
- Shannon
- Pielou evenness

Beta diversity:
- Bray-Curtis
- Jaccard

Faith's PD and UniFrac were not used for ITS.

## FUNGuild

Fixed commit:

```text
c2d1529394d77494d6c50e586df1e7b5b0a4e4df
```

Input:
- 5,272 fungal ASVs

Assigned:
- 4,204 ASVs

Unassigned:
- 1,068 ASVs

Confidence:
- Possible: 1,136
- Probable: 2,511
- Highly Probable: 557

Conservative subset:
- Probable + Highly Probable
- 3,068 ASVs
- 498,538 reads

Targeted ecological signatures:
- Any Saprotroph
- Endophyte
- Plant Pathogen
- Fungal Parasite
- Ericoid Mycorrhizal
- Ectomycorrhizal

Model:

```r
~ Fertilisation * PGPR + Cultivar + Bloc
```

## Interpretation

PICRUSt2 estimates predicted metagenomic functional potential and does not directly measure metabolic activity. FUNGuild provides taxonomically inferred ecological guilds and likewise does not directly measure biological activity.

Final biological interpretation should therefore distinguish:
- observed microbial community composition
- predicted/inferred functional potential
- directly measured biological activity
