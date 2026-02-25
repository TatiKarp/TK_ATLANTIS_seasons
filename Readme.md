# TK_ATLANTIS_seasons

Analysis of nasal brush RNAseq data from the **ATLANTIS** cohort study in relation to seasonal changes in gene expression, pathway activity, and immune cell composition.

---

## Overview

This repository contains the analysis pipeline for investigating how seasonal variation affects nasal brush transcriptomics in the ATLANTIS cohort. The analysis includes:

- Differential expression (DE) analysis using a 15-day seasonal window to define seasons
- DE analysis using Nov 15 cut off between 2 seasons
- Gene Set Variation Analysis (GSVA) of seasonal pathways
- Subgroup analyses in allergy and asthma
- Replication in the independent **PIAMA** cohort
- Viral transcript detection across seasons
- Cell type deconvolution across seasons
- Comparison with GWAS catalogue genes
- Pathway enrichment analysis
- Figure generation scripts

---

## Repository Structure

```
TK_ATLANTIS_seasons/
├── config.yml.example              # Template for user-specific paths and settings
│
├── 00_seasons_DE_15_shift.R        # DE analysis: seasonal shift (15-day window)
├── 01_season_15Nov_baseline_table.R # Baseline characteristics table (Nov 15 split)
├── 02_season_DE_15days.R           # DE analysis: 15-day seasonal windows
├── 03_GSVA_seasons_vs_date.R       # GSVA scores vs. collection date
├── 04a_seasonal_GSVA_in_allergy.R  # GSVA in allergy subgroup
├── 04b_seasonal_GSVA_in_asthma.R   # GSVA in asthma subgroup
├── 05a_PIAMA_DE_15Nov.R            # PIAMA replication: DE analysis
├── 05b_PIAMA_GSVA_season_15Nov.R   # PIAMA replication: GSVA
├── 06a_asthma_p_val_season_no_season.R  # Asthma seasonal vs. non-seasonal p-value comparison
├── 06b_compare genes_with_GWAS_catalogue.R  # Overlap with GWAS catalogue genes
├── 07a_viruses_in_seasons_15Nov.R  # Viral transcripts by season (R)
├── 07b_all_viruses_in_seasons_15Nov.py  # Viral transcripts by season (Python)
├── 08_season_15Nov_pathways.R      # Pathway enrichment across seasons
├── 09_cell_types_seasons_15Nov.R   # Cell type deconvolution by season
│
├── Figure_1_combined.R             # Figure 1 assembly
├── Figure_2_combined.R             # Figure 2 assembly
├── Figure_3_combined.R             # Figure 3 assembly
├── Figure_4_combined.R             # Figure 4 assembly
└── SF1_new_top_genes_vs_date.R     # Supplementary Figure 1: top genes vs. date
```

---

## Requirements

### R packages

The analysis requires R (≥ 4.2). Key packages include:

| Package | Purpose |
|--------|---------|
| `DESeq2` / `limma` | Differential expression analysis |
| `GSVA` | Gene set variation analysis |
| `edgeR` | Count normalisation |
| `ggplot2`, `patchwork` | Visualisation |
| `yaml` | Config file parsing |

> Install Bioconductor packages with:
> ```r
> if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
> BiocManager::install(c("DESeq2", "GSVA", "limma", "edgeR"))
> ```

### Python (script 07b)

```
python >= 3.8
pandas
numpy
matplotlib / seaborn
```

Install with:
```bash
pip install pandas numpy matplotlib seaborn
```

---

## Configuration

Copy the example config and fill in your local paths:

```bash
cp config.yml.example config.yml
```

Edit `config.yml` to set paths to your input data (count matrices, metadata, gene sets, etc.). **Do not commit `config.yml`** — it is already listed in `.gitignore`.

---

## Running the Analysis

Scripts are numbered to indicate the intended execution order. Run them sequentially from the repo root:

```r
source("00_seasons_DE_15_shift.R")
source("01_season_15Nov_baseline_table.R")
source("02_season_DE_15days.R")
# ... and so on
```

Figure scripts (`Figure_1_combined.R` etc.) depend on outputs from the numbered scripts and should be run last.

---

## Data Availability

The raw data from the ATLANTIS cohort are not included in this repository due to privacy and ethical restrictions. Access to the data can be requested via **[link or contact details]**.

The PIAMA cohort data used in scripts `05a` and `05b` are subject to a separate data access agreement.

---

## Citation

> **[Author(s)], [Year]. [Paper title]. [Journal]. DOI: [doi]**

---

## Contact

For questions about the analysis, please open a GitHub Issue or contact **[karp.tatiana.d@gmail.com/Tatiana Karp]**.

---
