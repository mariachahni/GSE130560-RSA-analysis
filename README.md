# GSE130560 Analysis: WP5637 Pathway Validation

This script was used to analyse the GSE130560 dataset as part of my bachelor's thesis at Maastricht University (BBS3006, 2026).

**Author:** Maria Chahni  
**Supervisor:** Dr. D. Slenter  
**Department:** Translational Genomics, FHML  

## Project Description

This repository contains the R script used to validate the computational pathway model WP5637 (WikiPathways), which models TGF-β1/SMAD2/3-mediated regulation of decidual NK cell cytotoxic activity. Validation was performed using publicly available single-cell RNA-seq data from GSE130560 (Wang et al. 2021), comparing gene expression in decidual NK cell subsets between RSA patients and healthy controls.

## Main Steps

1. Load and prepare GSE130560 into a Seurat object
2. Log-normalize gene expression data
3. Run differential expression analysis per dNK subset (dNKp, dNK1-4) using Wilcoxon rank sum test
4. Perform over-representation analysis (ORA) using Fisher's exact test against WP5637 pathway nodes
5. Generate heatmap of pathway node expression across dNK subsets
6. Run WikiPathways-wide ORA using clusterProfiler with Benjamini-Hochberg correction
7. Export dNK4 expression data for PathVisio overlay

## How to Run

1. Download the dataset from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130560
   - GSE130560_phenotype.csv.gz
   - GSE130560_matrix.RData.gz
2. Place both files in your working directory
3. Update the `setwd()` path in Section 2 to your working directory
4. Run the full script top to bottom in R (version 4.6.0)
5. Install required packages if not already available:
   - Seurat, ggplot2, ggrepel, R.utils, clusterProfiler, rWikiPathways, org.Hs.eg.db, enrichplot, dplyr, tidyr

## Output Files

- `ORA_summary.csv` : Fisher's exact test results per dNK subset
- `heatmap_pathway_nodes_final.png` : Pathway node expression heatmap
- `WikiPathways_ORA_[subset].csv` : WikiPathways ORA results per subset
- `WikiPathways_ORA_dotplot.png` : ORA dotplot
- `pathvisio_expression_dNK4.txt` : Expression data for PathVisio overlay


## Related Resources

- Pathway WP5637: https://www.wikipathways.org/pathways/WP5637.html
- Dataset: Wang et al. 2021, https://doi.org/10.1016/j.gpb.2020.11.002
