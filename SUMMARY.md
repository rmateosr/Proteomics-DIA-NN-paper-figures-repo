# Summary of R Scripts for DIA-NN Proteomics Paper Figures

**Repository:** https://github.com/rmateosr/Proteomics-DIA-NN-paper-figures-repo
**Date:** 2026-02-10
**Total scripts:** 17

---

## What was done

1. Reviewed all 17 R scripts located in `C:\Users\Raul\Dropbox\Papers\DIA-NN\Figures\scripts`
2. Created a public GitHub repository at https://github.com/rmateosr/Proteomics-DIA-NN-paper-figures-repo
3. Initialized the working scripts directory as a git repository linked to the remote GitHub repo
4. Committed and pushed all 17 scripts to the `main` branch

---

## Script Categories

The scripts fall into five groups:

### 1. Hotspot Mutation Visualization (Cell Lines)
These scripts detect and visualize cancer hotspot mutations at the peptide level in a panel of 10 cell lines (HeLa, A549, CCRF-CEM, COLO205, H226, H23, RPMI8226, T47D, 4pool, NCI7ref) analyzed in duplicate.

| Script | Description |
|--------|-------------|
| `Hotspot_figure_script_cellline.R` | Exploratory figures for KRAS, TP53, and BRAF hotspot mutations |
| `Hotspot_figure_script_cellline_forpaper.R` | Paper-quality versions (KRAS, TP53, BRAF, TP53-T47D) with refined formatting |
| `Hotspot_figure_script_cellline_T47D.R` | Focused on TP53 L194F homozygous mutation in T47D cell line |
| `Hotspot_KEAP1_WT_figure_script_cellline_forWakate.R` | KEAP1 wild-type protein signal across cell lines (Wakate presentation) |
| `Hotspot_KRAS_G12A_figure_script_cellline_forWakate.R` | KRAS G12A mutation peptide detection (Wakate presentation) |
| `Hotspot_KRAS_G12C_figure_script_cellline_forWakate.R` | KRAS G12C mutation peptide detection (Wakate presentation) |

**Key method:** All use a canonical peptide re-integration algorithm that finds the wild-type counterpart of each mutated peptide, handling edge cases where mutations create or destroy tryptic cleavage sites (R/K residues).

### 2. Hotspot Mutation Visualization (PDX Samples)
| Script | Description |
|--------|-------------|
| `Hotspot_figure_script_PDX_forpaper.R` | Visualizes KRAS G12/G13, TP53, and CTNNB1 mutations across 49 PDX tumor samples |

### 3. Gene Fusion Analysis
| Script | Description |
|--------|-------------|
| `noncanonicalpeptidesanalysis_GeneFusion_forpaper.R` | Comprehensive gene fusion peptide analysis across all detected fusions |
| `noncanonicalpeptidesanalysis_GeneFusion_JustGAS6_forpaper.R` | Focused on GAS6-RASA3 gene fusion peptide |
| `SummaryofGeneFusionsbyFusionPDB_Adachisamples.R` | OncoPrint visualizations of gene fusions at three evidence levels (Level 1, 2, 3) |
| `SummaryofGeneFusionsbyFusionPDB_Adachisamples_forpaper.R` | Paper-quality OncoPrint (Level 1 only, 1000 DPI) |

### 4. Protein Expression
| Script | Description |
|--------|-------------|
| `PREX2_PDX_Geneexpression_forpaper.R` | PREX2 protein signal across PDX samples (paper version) |
| `PREX2_protein_expression_figure_script.R` | PREX2 protein expression in PDX (larger figure dimensions) |

### 5. Hotspot Mutation Summary & NRF2/KEAP1 Classifier
| Script | Description |
|--------|-------------|
| `Summaryofhotspotmutations_Adachisamples.R` | OncoPrint of hotspot mutations from WES data across PDX samples |
| `Summaryofhotspotmutations_Adachisamples_forpaper.R` | Paper-quality version of the hotspot mutation OncoPrint |
| `SJpaperfigurereference.R` | NRF2/KEAP1 splicing junction classifier: boxplots, PR-AUC curves, cross-validation (most complex script, ~925 lines) |
| `SJpaperfigurereference - Copy.R` | Copy of the above |

---

## Common Patterns Across Scripts

- **Normalization:** Column-sum normalization (value / column_sum x 1,000,000) used throughout
- **Data granularities:** `pr_matrix` (precursor level, for mutation-specific peptides) and `pg_matrix` (protein group level, for protein abundance)
- **Sample types:** Cell line panel (10 lines in duplicate) and PDX cohort (49 tumors)
- **Versioning:** Many scripts exist as exploratory and `_forpaper` pairs with refined formatting
- **Key libraries:** ggplot2, dplyr, data.table, stringr, reshape2, ComplexHeatmap, PRROC, cowplot

---

## Notes

- `SJpaperfigurereference - Copy.R` is an identical duplicate of `SJpaperfigurereference.R`
- `Summaryofhotspotmutations_Adachisamples_forpaper.R` contains an incomplete `write.table()` call on line 127 that would cause an error if executed
