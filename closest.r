#install.packages("DescTools")
library(DescTools)
tss_gene = c(myProm$tssStart)
tss_transcript = c(gene_data$transcript_start_site)
probe_start = c(myProm$probeStart)
closest_tss = list()

c = lapply(probe_start, function(i) (Closest(probe_start, tss_gene)))

closest_probe = unique(unlist(c))

#closest_t = lapply(tss_transcript, function(i) (Closest(tss_transcript, tss_gene)))
#closest_transcript = unique(unlist(closest_t))
#test = gene_data[gene_data$transcript_start_site == "8076939",]


myProm = as.data.frame(myProm)
test = myProm[myProm$probeStart == closest_probe, , ]
probe = test$probes

probe_colors <- rep(c("lightblue", "lightgreen"), each = 5)

ra = rowAnnotation(probe = anno_block(labels = probe_start,gp = gpar(fill = 2), labels_gp = gpar(fontsize = 2)))
ra =  rowAnnotation(foo = anno_block(
         panel_fun = function(index, levels) {
                 grid.rect(gp = gpar(fill = col, col = "black"))
                 grid.text(paste(levels, collapse = ","), 0.5, 0.5, rot = 90,
                         gp = gpar(col = col[levels[1]]))}))

h2 = draw(Heatmap(beta.data,
             col = my_palette, 
             name ="Methylation degree",
             cluster_rows = FALSE,
             column_order = column_order(h1),
             clustering_method_columns = "ward.D2",
             column_title = paste(my.gene,"promoter: upS=",my.upS,"dnS=",my.dnS, "\n\n\n"),
             show_heatmap_legend = TRUE,
             column_title_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 0.5),
             row_names_gp = gpar(fontsize = 5),
             width = 5,
             left_annotation = ra))
