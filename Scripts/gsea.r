###

# Description - performs gene enerichment analysis

# Load necessary libraries
library(msigdbr)
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
load("variance_percent.Rdata") # from filtering.R script

gene_data = vd[,-c(1:2,4,6:length(vd))] %>%
  distinct(gene_symbol, .keep_all = TRUE)  # Remove duplicates # remove other attributes except names needed to perform gsea

# Get hallmark gene sets for Homo sapiens
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(grepl("HALLMARK_", gs_name)) %>%
  dplyr::select(gs_name, gene_symbol)

# Create named vector for gene list and sort
gene_list_named <- setNames(gene_data$transcript_variance, gene_data$gene_symbol)
gene_list_named <- sort(gene_list_named, decreasing = TRUE)

# Perform GSEA with hallmark gene sets
gsea_results <- GSEA(
  geneList = gene_list_named,
  TERM2GENE = hallmark_sets,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.4,
  pAdjustMethod = "BH",
  scoreType = "pos"  # Use positive score type
)

# Plot results
dotplot(gsea_results, showCategory = 20) + 
  ggtitle("GSEA Dot Plot for Hallmark Gene Sets") + 
  theme_minimal()

# for other pathways
gene_symbols = unique(gene_data$gene_symbol)
gene_ids = bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_enrich = enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
# change the ontology option to MF for molecular function, BP for Biological function and CC for Cellular component

plot_go = dotplot(go_enrich, showCategory = 15)
barplot(go_enrich, showCategory = 15)
ggsave("go_enrich_plot.png", plot = plot_go, width = 10, height = 6)
