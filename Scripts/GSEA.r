###

# Description - performs gene enerichment analysis

# Load necessary libraries
library(msigdbr)
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
load("variance_percent.Rdata") # from filtering.R script

gene_data = vd[,-c(1:2,4,6:length(vd))] # remove other attributes except names needed to perform gsea

gene_symbols = unique(gene_data$gene_symbol)
gene_ids = bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_enrich = enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
# change the ontology option to MF for molecular function, BP for Biological function and CC for Cellular component

plot_go = dotplot(go_enrich, showCategory = 15)
barplot(go_enrich, showCategory = 15)
ggsave("go_enrich_plot.png", plot = plot_go, width = 10, height = 6)

