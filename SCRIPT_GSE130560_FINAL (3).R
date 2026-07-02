# =============================================================
# GSE130560 Complete Analysis Script
# TGF-β1/SMAD2/3 pathway (WP5637) in decidual NK cells: RSA vs Healthy
# Dataset: Wang et al. 2021, Genomics Proteomics & Bioinformatics
# Pathway: WP5637, WikiPathways
# Analysis by: Maria Chahni
# =============================================================

# ---- 1. Load required libraries ----
library(Seurat)
library(ggplot2)
library(ggrepel)
library(R.utils)
library(clusterProfiler)
library(rWikiPathways)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(tidyr)

# ---- 2. Set working directory ----
# NOTE: update this path to the folder containing the downloaded
# GSE130560 files (GSE130560_phenotype.csv.gz and GSE130560_matrix.RData.gz)
setwd("C:/Users/gaber/Downloads")

# ---- 3. Load dataset ----
# Dataset: GSE130560 (Wang et al. 2021, Genomics Proteomics Bioinformatics)
# Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130560
# If files are still gzipped (first run), unzip first by uncommenting:
#   GSE130560 <- read.csv(gunzip("GSE130560_phenotype.csv.gz"))
#   load(gunzip("GSE130560_matrix.RData.gz"))
# If already unzipped (subsequent runs), load directly:
GSE130560 <- read.csv("GSE130560_phenotype.csv")
load("GSE130560_matrix.RData")

# ---- 4. Create and prepare Seurat object ----
GSE130560_SeuratObject <- CreateSeuratObject(counts = matrix)
rownames(GSE130560) <- GSE130560$X
GSE130560_SeuratObject <- AddMetaData(GSE130560_SeuratObject, metadata = GSE130560)
GSE130560_SeuratObject <- NormalizeData(GSE130560_SeuratObject)

# Check cell types and disease groups
table(GSE130560_SeuratObject$celltype)
table(GSE130560_SeuratObject$disease)

# ---- 5. Define WP5637 pathway nodes ----
# Gene symbols standardised to HGNC nomenclature
# (e.g. IL8 -> CXCL8, NKG2D -> KLRK1)
pathway_nodes <- c("CD9", "SMAD4", "ANGPT1", "TGFB1", "SMAD3", "KIR2DL4",
                    "IFNG", "SMURF2", "CD103", "NKG2D", "CD14", "TGFBR2",
                    "TGFBR1", "SMAD7", "ANGPT2", "VEGFC", "PLGF", "SMAD2",
                    "CXCL8", "CXCL10", "KLRK1", "FCGR3A", "IL2RA", "FOXP3", "PGF")

# ---- 6. Differential expression analysis per dNK subset ----
# Wilcoxon rank sum test (Seurat::FindMarkers default), RSA vs CTRL

# dNKp (dNK_a)
dNKa <- subset(GSE130560_SeuratObject, subset = celltype == "dNK_a")
Idents(dNKa) <- "disease"
markers_dNKa <- FindMarkers(dNKa, ident.1 = "RSA", ident.2 = "CTRL")
markers_dNKa$gene <- rownames(markers_dNKa)
markers_dNKa$subset <- "dNKp"

# dNK1 (dNK_b)
dNKb <- subset(GSE130560_SeuratObject, subset = celltype == "dNK_b")
Idents(dNKb) <- "disease"
markers_dNKb <- FindMarkers(dNKb, ident.1 = "RSA", ident.2 = "CTRL")
markers_dNKb$gene <- rownames(markers_dNKb)
markers_dNKb$subset <- "dNK1"

# dNK2 (dNK_c)
dNKc <- subset(GSE130560_SeuratObject, subset = celltype == "dNK_c")
Idents(dNKc) <- "disease"
markers_dNKc <- FindMarkers(dNKc, ident.1 = "RSA", ident.2 = "CTRL")
markers_dNKc$gene <- rownames(markers_dNKc)
markers_dNKc$subset <- "dNK2"

# dNK3 (dNK_d)
dNKd <- subset(GSE130560_SeuratObject, subset = celltype == "dNK_d")
Idents(dNKd) <- "disease"
markers_dNKd <- FindMarkers(dNKd, ident.1 = "RSA", ident.2 = "CTRL")
markers_dNKd$gene <- rownames(markers_dNKd)
markers_dNKd$subset <- "dNK3"

# dNK4 (dNK_e)
dNKe <- subset(GSE130560_SeuratObject, subset = celltype == "dNK_e")
Idents(dNKe) <- "disease"
markers_dNKe <- FindMarkers(dNKe, ident.1 = "RSA", ident.2 = "CTRL")
markers_dNKe$gene <- rownames(markers_dNKe)
markers_dNKe$subset <- "dNK4"

# ---- 7. Over-Representation Analysis (ORA): Fisher's exact test ----
# Tests whether WP5637 pathway nodes are over-represented among
# significant DEGs (Bonferroni-adjusted p < 0.05) per dNK subset

run_ORA <- function(markers, pathway_nodes) {
  overlap <- pathway_nodes[pathway_nodes %in% rownames(markers)]
  sig_degs <- rownames(markers[markers$p_val_adj < 0.05, ])
  total_genes <- nrow(markers)
  pathway_in_sig <- sum(overlap %in% sig_degs)
  pathway_not_sig <- length(overlap) - pathway_in_sig
  non_pathway_sig <- length(sig_degs) - pathway_in_sig
  non_pathway_not_sig <- total_genes - pathway_not_sig - non_pathway_sig - pathway_in_sig
  contingency_table <- matrix(c(pathway_in_sig, pathway_not_sig,
                                 non_pathway_sig, non_pathway_not_sig), nrow = 2)
  fisher.test(contingency_table, alternative = "greater")
}

fisher_a <- run_ORA(markers_dNKa, pathway_nodes)  # dNKp
fisher_b <- run_ORA(markers_dNKb, pathway_nodes)  # dNK1
fisher_c <- run_ORA(markers_dNKc, pathway_nodes)  # dNK2
fisher_d <- run_ORA(markers_dNKd, pathway_nodes)  # dNK3
fisher_e <- run_ORA(markers_dNKe, pathway_nodes)  # dNK4

# View results
fisher_a; fisher_b; fisher_c; fisher_d; fisher_e

# ---- 8. ORA summary table ----
ora_summary <- data.frame(
  Subset = c("dNKp", "dNK1", "dNK2", "dNK3", "dNK4"),
  Corresponds_to = c("dNK_a", "dNK_b", "dNK_c", "dNK_d", "dNK_e"),
  Total_DEGs = c(
    nrow(markers_dNKa[markers_dNKa$p_val_adj < 0.05,]),
    nrow(markers_dNKb[markers_dNKb$p_val_adj < 0.05,]),
    nrow(markers_dNKc[markers_dNKc$p_val_adj < 0.05,]),
    nrow(markers_dNKd[markers_dNKd$p_val_adj < 0.05,]),
    nrow(markers_dNKe[markers_dNKe$p_val_adj < 0.05,])
  ),
  p_value = round(c(fisher_a$p.value, fisher_b$p.value, fisher_c$p.value,
                     fisher_d$p.value, fisher_e$p.value), 4),
  Odds_Ratio = round(c(fisher_a$estimate, fisher_b$estimate, fisher_c$estimate,
                        fisher_d$estimate, fisher_e$estimate), 2)
)
write.csv(ora_summary, "ORA_summary.csv", row.names = FALSE)

# ---- 9. Heatmap of pathway node expression across dNK subsets ----

get_val <- function(markers, genes, col) {
  sapply(genes, function(g) {
    if (g %in% rownames(markers)) markers[g, col] else NA
  })
}

overlap_all <- unique(c(
  pathway_nodes[pathway_nodes %in% rownames(markers_dNKa)],
  pathway_nodes[pathway_nodes %in% rownames(markers_dNKb)],
  pathway_nodes[pathway_nodes %in% rownames(markers_dNKc)],
  pathway_nodes[pathway_nodes %in% rownames(markers_dNKd)],
  pathway_nodes[pathway_nodes %in% rownames(markers_dNKe)]
))

heatmap_data <- data.frame(
  gene = rep(overlap_all, 5),
  subset = rep(c("dNKp", "dNK1", "dNK2", "dNK3", "dNK4"), each = length(overlap_all)),
  log2FC = c(get_val(markers_dNKa, overlap_all, "avg_log2FC"),
             get_val(markers_dNKb, overlap_all, "avg_log2FC"),
             get_val(markers_dNKc, overlap_all, "avg_log2FC"),
             get_val(markers_dNKd, overlap_all, "avg_log2FC"),
             get_val(markers_dNKe, overlap_all, "avg_log2FC")),
  significant = c(
    ifelse(overlap_all %in% rownames(markers_dNKa[markers_dNKa$p_val_adj < 0.05,]), "*", ""),
    ifelse(overlap_all %in% rownames(markers_dNKb[markers_dNKb$p_val_adj < 0.05,]), "*", ""),
    ifelse(overlap_all %in% rownames(markers_dNKc[markers_dNKc$p_val_adj < 0.05,]), "*", ""),
    ifelse(overlap_all %in% rownames(markers_dNKd[markers_dNKd$p_val_adj < 0.05,]), "*", ""),
    ifelse(overlap_all %in% rownames(markers_dNKe[markers_dNKe$p_val_adj < 0.05,]), "*", "")
  )
)

heatmap_data$subset <- factor(heatmap_data$subset,
                               levels = c("dNKp", "dNK1", "dNK2", "dNK3", "dNK4"))

# Heatmap with numeric log2FC labels and bold border for significant tiles
heatmap_data$log2FC_label <- ifelse(is.na(heatmap_data$log2FC), "", round(heatmap_data$log2FC, 2))
heatmap_data$border_width <- ifelse(heatmap_data$significant == "*", 1.5, 0.3)

ggplot(heatmap_data, aes(x = subset, y = gene, fill = log2FC)) +
  geom_tile(aes(color = significant == "*"), linewidth = heatmap_data$border_width) +
  geom_text(aes(label = log2FC_label), size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "log2FC", na.value = "grey90") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  labs(title = "Pathway node expression across dNK subsets: RSA vs Healthy",
       x = "dNK Subset", y = "Pathway Gene",
       caption = "Bold black border indicates Bonferroni-adjusted p < 0.05 | Grey = not detected") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))

ggsave("heatmap_pathway_nodes_final.png", width = 10, height = 8, dpi = 300)

# ---- 10. WikiPathways Over-Representation Analysis ----
# Tests whether DEGs per subset are enriched in any WikiPathways
# pathway (database-wide), using Benjamini-Hochberg correction

wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens", format = "gmt")
wp2gene <- clusterProfiler::read.gmt(wp.hs.gmt)
wp2gene <- wp2gene %>% tidyr::separate(term, c("name", "version", "wpid", "org"), sep = "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)
wpid2name <- wp2gene %>% dplyr::select(wpid, name)

# Convert significant DEGs (gene symbols) to Entrez IDs
get_entrez_degs <- function(markers) {
  sig_genes <- rownames(markers[markers$p_val_adj < 0.05, ])
  if (length(sig_genes) == 0) return(NULL)
  entrez <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(entrez$ENTREZID)
}

# Background gene universe (all genes tested in dNK4; used for all subsets)
bg_genes <- bitr(rownames(markers_dNKe), fromType = "SYMBOL",
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

degs_a <- get_entrez_degs(markers_dNKa)
degs_b <- get_entrez_degs(markers_dNKb)
degs_c <- get_entrez_degs(markers_dNKc)
degs_d <- get_entrez_degs(markers_dNKd)
degs_e <- get_entrez_degs(markers_dNKe)

run_WP_ORA <- function(degs, bg, name) {
  if (is.null(degs) || length(degs) < 2) {
    message(paste("Skipping", name, "- too few DEGs"))
    return(NULL)
  }
  enrichWP(gene = degs, universe = bg, organism = "Homo sapiens",
           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
}

ora_dNKp <- run_WP_ORA(degs_a, bg_genes, "dNKp")
ora_dNK1 <- run_WP_ORA(degs_b, bg_genes, "dNK1")
ora_dNK2 <- run_WP_ORA(degs_c, bg_genes, "dNK2")
ora_dNK3 <- run_WP_ORA(degs_d, bg_genes, "dNK3")
ora_dNK4 <- run_WP_ORA(degs_e, bg_genes, "dNK4")

# Save significant WikiPathways ORA results per subset
save_ora_results <- function(ora, name) {
  if (!is.null(ora) && nrow(ora@result) > 0) {
    sig <- ora@result[ora@result$p.adjust < 0.05, ]
    if (nrow(sig) > 0) {
      write.csv(sig, paste0("WikiPathways_ORA_", name, ".csv"), row.names = FALSE)
    }
  }
}
save_ora_results(ora_dNKp, "dNKp")
save_ora_results(ora_dNK1, "dNK1")
save_ora_results(ora_dNK2, "dNK2")
save_ora_results(ora_dNK3, "dNK3")
save_ora_results(ora_dNK4, "dNK4")

# ---- 11. WikiPathways ORA dotplot ----
combine_ora <- function(ora, subset_name) {
  if (is.null(ora)) return(NULL)
  sig <- ora@result[ora@result$p.adjust < 0.05, ]
  if (nrow(sig) == 0) return(NULL)
  sig$subset <- subset_name
  return(sig)
}

all_ora <- rbind(
  combine_ora(ora_dNKp, "dNKp"), combine_ora(ora_dNK1, "dNK1"),
  combine_ora(ora_dNK2, "dNK2"), combine_ora(ora_dNK3, "dNK3"),
  combine_ora(ora_dNK4, "dNK4")
)

# Restrict to biologically relevant immune/vascular pathways for the main figure
keep_pathways <- c(
  "Allograft rejection",
  "Inflammatory response pathway",
  "Uterine natural killer cells and progesterone estrogen and chorionic gonadotropin",
  "Th17 cell differentiation pathway",
  "Prostaglandin signaling",
  "Overview of proinflammatory and profibrotic mediators",
  "VEGFA VEGFR2 signaling",
  "Nonalcoholic fatty liver disease",
  "Electron transport chain OXPHOS system in mitochondria",
  "Modulators of TCR signaling and T cell activation"
)

filtered_ora <- all_ora[all_ora$Description %in% keep_pathways, ]
filtered_ora$Description <- gsub(
  "Uterine natural killer cells and progesterone estrogen and chorionic gonadotropin",
  "Uterine NK cells and hormones (WP5569)", filtered_ora$Description)
filtered_ora$subset <- factor(filtered_ora$subset,
                               levels = c("dNKp", "dNK1", "dNK2", "dNK3", "dNK4"))

p <- ggplot(filtered_ora, aes(x = subset, y = Description, size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#E85D4A", high = "#A8C5DA", name = "Adjusted\np-value") +
  scale_size_continuous(range = c(4, 12), name = "Gene Count") +
  labs(title = "WikiPathways ORA across dNK subsets: RSA vs Healthy",
       x = "dNK Subset", y = "") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 9),
        plot.title = element_text(size = 11), plot.margin = margin(10, 20, 10, 10))

ggsave("WikiPathways_ORA_dotplot.png", plot = p, width = 12, height = 7, dpi = 300)

# ---- 12. PathVisio expression overlay export (dNK4) ----
# Exports log2FC and significance for WP5637 nodes detected in dNK4,
# for import into PathVisio (Data > Select expression data) to
# visualise expression directly on the pathway model.

overlap <- pathway_nodes[pathway_nodes %in% rownames(markers_dNKe)]
pathvisio_data <- data.frame(
  GeneID = overlap,
  log2FC = round(markers_dNKe[overlap, "avg_log2FC"], 3),
  p_adj = markers_dNKe[overlap, "p_val_adj"],
  significant = ifelse(markers_dNKe[overlap, "p_val_adj"] < 0.05, "Yes", "No")
)
write.table(pathvisio_data, "pathvisio_expression_dNK4.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# =============================================================
# End of script
# Outputs generated:
#   ORA_summary.csv
#   heatmap_pathway_nodes_final.png
#   WikiPathways_ORA_dNKp.csv ... WikiPathways_ORA_dNK4.csv
#   WikiPathways_ORA_dotplot.png
#   pathvisio_expression_dNK4.txt
# =============================================================
