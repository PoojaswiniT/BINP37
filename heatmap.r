h1 = (Heatmap(gene_matrix,
        col = cols_exp, 
        name ="Expression levels", 
        column_title = paste0("Expression of ", unique(gene_data$gene_symbol), " gene ", "strand: ", strand),
        cluster_columns = TRUE,
        row_dend_width = unit(0.5, "cm"),
        show_heatmap_legend = TRUE,
        clustering_method_columns = "ward.D2",
        column_title_gp = gpar(fontsize = 10),
        column_names_gp =  gpar(fontsize = 0.5),
        row_names_gp = gpar(fontsize = 5),
        #top_annotation = column_annotation,
        width = 12,
        column_km = 4))

col_order_h1 <- as.integer(unlist(column_order(h1)))
########### No TSS Annotation########
h2 <- (Heatmap(beta.data,
                   col = my_palette, 
                   name ="Methylation degree",
                   cluster_rows = FALSE,
                   column_order = col_order_h1,
                   clustering_method_columns = "ward.D2",
                   column_title = paste(my.gene,"promoter: upS=",my.upS,"dnS=",my.dnS, "\n"),
                   show_heatmap_legend = TRUE,
                   column_title_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 0.5),
                   row_names_gp = gpar(fontsize = 5),
                   width = 12))
col_order_h2 <- column_order(h2)

# Check if column orders are identical 
stopifnot(identical((col_order_h1), (col_order_h2)))

# Plot both heatmaps together
#library(ggplot2)
library(gridExtra)
#library(plotly)
# Arrange the heatmaps in a grid
grob_h1 = ggplotify::as.grob(h1)
grob_h2 = ggplotify::as.grob(h2)

# Arrange the grobs in a grid
grid.arrange(grob_h1, grob_h2, ncol = 1)
######## With TSS Annotation acc to gene#######
gd = cg_data_arranged[cg_data_arranged$gene_symbol == "ATP6V0E1", ]

tss_gene = c(unique(gd$tssStart))
probe_start = c(unique(gd$probeStart))
enst_ids = c(unique(gd$transcript_id))
enst_tss = c(unique(gd$transcript_start_site))

library(ComplexHeatmap)
library(grid)
#pdf("heatmap_output.pdf")

# Loop through each tss_gene value
for (tss_value in tss_gene) {
  # Find the index of the closest probeStart to the current tss_gene value
  closest_probe_index <- which.min(abs(probe_start - tss_value))
  closest_tss_index <- which.min(abs(enst_tss - tss_value))
  # Get the closest probeStart value using the index
  closest_probe_start <- probe_start[closest_probe_index]
  closest_enst_tss <- enst_tss[closest_tss_index]
  probe_id <- unique(gd$probes[gd$probeStart == closest_probe_start])
  enst_id <- unique(gd$transcript_id[gd$transcript_start_site == closest_enst_tss])
  # Find the column index (probe_pos) for the probe in the heatmap
  probe_pos <- which(rownames(beta.data) == probe_id)
  enst_pos <- which(rownames(gene_matrix) == enst_id)
  # Annotate heatmap with color change
  probe_col <- rep("white", nrow(beta.data))
  probe_col[probe_pos] <- "black"
}

# Create a named vector for colors
color_vector <- c("black", "white")
names(color_vector) <- c("black", "white")

h2 <- (Heatmap(beta.data,
              col = my_palette, 
              name = "Methylation degree",
              cluster_rows = FALSE,
              column_order = col_order_h1,
              clustering_method_columns = "ward.D2",
              column_title = paste(my.gene, "promoter: upS=", my.upS, "dnS=", my.dnS, "\n\n\n"),
              show_heatmap_legend = TRUE,
              column_title_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 0.5),
              row_names_gp = gpar(fontsize = 5),
              width = 12,
              row_title_side = "left",  # Place row titles on the left side
              left_annotation = HeatmapAnnotation(df = data.frame(probe_col), 
                                                  col = list(probe_col = color_vector),
                                                  which = "row", show_annotation_name = FALSE, show_legend = FALSE)))
col_order_h2 <- column_order(h2)

# Check if column orders are identical 
stopifnot(identical((col_order_h1), (col_order_h2)))

# Plot both heatmaps together
#library(ggplot2)
library(gridExtra)
#library(plotly)
# Arrange the heatmaps in a grid
grob_h1 = ggplotify::as.grob(h1)
grob_h2 = ggplotify::as.grob(h2)

# Arrange the grobs in a grid
grid.arrange(grob_h1, grob_h2, ncol = 1)
######## TSS annotation acc to ENSTs #######


gd = cg_data_arranged[cg_data_arranged$gene_symbol == "ESR1", ]

tss_gene = c(unique(gd$tssStart))
probe_start = c(unique(gd$probeStart))
enst_ids = c(unique(gd$transcript_id))
enst_tss = c(unique(gd$transcript_start_site))

library(ComplexHeatmap)
library(grid)
probe_pos <- numeric(length(enst_tss))
closest_probe_ids <- character(length(enst_tss))
probe_col <- rep("white", nrow(beta.data))

for (i in seq_along(enst_tss)) {
  tss_value <- enst_tss[i]
  # Find the index of the closest probeStart to the current tss_gene value
  closest_probe_index <- which.min(abs(probe_start - tss_value))
  # Get the closest probeStart value using the index
  closest_probe_start <- probe_start[closest_probe_index]
  probe_id <- unique(gd$probes[gd$probeStart == closest_probe_start])
  probe_pos[i] <- which(rownames(beta.data) == probe_id)
  closest_probe_ids[i] <- probe_id
  probe_col[probe_pos[i]] <- "black"
}

color_vector <- c("black", "white")
names(color_vector) <- c("black", "white")

h2 <- (Heatmap(beta.data,
              col = my_palette, 
              name = "Methylation degree",
              cluster_rows = FALSE,
              column_order = col_order_h1,
              clustering_method_columns = "ward.D2",
              column_title = paste(my.gene, "promoter: upS=", my.upS, "dnS=", my.dnS, "\n\n\n"),
              show_heatmap_legend = TRUE,
              column_title_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 0.5),
              row_names_gp = gpar(fontsize = 5),
              width = 12,
              row_title_side = "left",  # Place row titles on the left side
              left_annotation = HeatmapAnnotation(df = data.frame(probe_col), 
                                                  col = list(probe_col = color_vector),
                                                  which = "row", show_annotation_name = FALSE, show_legend = FALSE)))
col_order_h2 <- column_order(h2)

# Check if column orders are identical 
stopifnot(identical((col_order_h1), (col_order_h2)))

# Plot both heatmaps together
#library(ggplot2)
library(gridExtra)
#library(plotly)
# Arrange the heatmaps in a grid
grob_h1 = ggplotify::as.grob(h1)
grob_h2 = ggplotify::as.grob(h2)

# Arrange the grobs in a grid
grid.arrange(grob_h1, grob_h2, ncol = 1)
