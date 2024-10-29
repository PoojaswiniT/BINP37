###

#Description - plotting box plots to visualise the correlation between subtypes and transcript expression


##### Load libraries #####

library(ggplot2)
library(reshape2)
library(gridExtra)

#Load the data
load("variance_data.Rdata") # From filtering.R script
for( i in unique(vd$gene_symbol)){
  gene_data = vd[vd$gene_symbol == i, ]
  # Remove the first two columns
  gene_count_only = gene_data[,-c(1:12)]

  gene_matrix = as.matrix(gene_count_only)
  colnames(gene_matrix) = sapply(sample_name_split(colnames(gene_matrix)), `[`, 1)
  rownames(gene_matrix) = c(gene_data$transcript_id)

  gene_matrix_transposed = t(gene_matrix)


  # Combine gene_matrix_transposed with NCN_levels
  combined_data = cbind(gene_matrix_transposed, NCN = NCN_levels)

  # Convert to data frame
  combined_data_df = as.data.frame(combined_data)

  # Melt the data for plotting
  melted_data = melt(combined_data_df, id.vars = "NCN", variable.name = "ENST_ID", value.name = "Count")

  # List to store plots
  plots_list = list()

  #colors_list = list()
  col_list = c("1"= "red", "2"= "purple", "3"= "blue", "4" = "lightblue", "5"= "green","6"= "gray")

  # Loop through each ENST_ID and create a plot
  for (enst in unique(melted_data$ENST_ID)) {
    # Subset data for the current ENST_ID
    enst_data = melted_data[melted_data$ENST_ID== enst, ]
    # Create plot for the current ENST_ID
    p = ggplot(enst_data, aes(x = factor(NCN_levels), y = Count, fill = factor(NCN)), xlim(c(0.5,0.3)), ylim(0,35)) +
      geom_boxplot() +
      labs(x = "Tumor subtype", y = "Count", title = paste("ENST ID:", enst, "  ")) +
      scale_fill_manual(values = col_list) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Store the plot in the list
    plots_list[[enst]] = p
  }

  #pdf(paste("D:/BINP37/box_plots/", "i".pdf")) # Change the path appropriately

  # Combine all plots into one
  combined_plot = do.call(grid.arrange, c(plots_list, ncol = 3))
  dev.off()
}
