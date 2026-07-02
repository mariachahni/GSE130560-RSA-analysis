# GSE130560-RSA-analysis
R analysis pipeline for thesis: TGF-β/SMAD2/3 regulation of dNK cytotoxicity
# GSE130560 Analysis — WP5637 Pathway Validation

This repository contains the R analysis script used for my bachelor's thesis:
**"TGF-β/SMAD2/3-Mediated Regulation of Decidual Natural Killer Cell Cytotoxicity 
at the Maternal-Fetal Interface and Its Disruption in Recurrent Pregnancy Loss"**

Author: Maria Chahni | Maastricht University | BBS3006 Internship 2026  
Supervisor: Dr. D. Slenter

---

## Project Description
This script validates the computational pathway model WP5637 (WikiPathways) 
using single-cell RNA-seq data from GSE130560 (Wang et al. 2021), comparing 
decidual NK cell gene expression between RSA patients and healthy controls.

---

## Main Steps
1. Load and prepare the GSE130560 dataset into a Seurat object
2. Log-normalize gene expression data
3. Run differential expression analysis (Wilcoxon rank sum test) per dNK subset (dNKp, dNK1-4)
4. Perform over-representation analysis (ORA) using Fisher's exact test against WP5637 pathway nodes
5. Generate heatmap of pathway node expression across dNK subsets
6. Run WikiPathways-wide ORA using clusterProfiler with Benjamini-Hochberg correction
7. Export expression data for PathVisio overlay (dNK4 subset)

---

## How to Run
1. Download the dataset from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130560
   - GSE130560_phenotype.csv.gz
   - GSE130560_matrix.RData.gz
2. Place both files in your working directory
3. Update the `setwd()` path in Section 2 to match your working directory
4. Run the script from top to bottom in R (version 4.6.0)
5. Required packages: Seurat, ggplot2, ggrepel, R.utils, clusterProfiler, 
   rWikiPathways, org.Hs.eg.db, enrichplot, dplyr, tidyr

---

## Output Files
- `ORA_summary.csv` — Fisher's exact test results per dNK subset
- `heatmap_pathway_nodes_final.png` — Pathway node expression heatmap
- `WikiPathways_ORA_[subset].csv` — WikiPathways ORA results per subset
- `WikiPathways_ORA_dotplot.png` — ORA dotplot figure
- `pathvisio_expression_dNK4.txt` — Expression overlay file for PathVisio

---

## Related Resources
- Pathway WP5637: https://www.wikipathways.org/pathways/WP5637.html
- Dataset paper: Wang et al. 2021, https://doi.org/10.1016/j.gpb.2020.11.002
