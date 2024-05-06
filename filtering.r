library(dplyr)
library(tidyr)

# Count number of transcripts
# Group by gene id, count the number of distinct transcript ids for each gene
transcript_count_per_gene = count_mapped_data %>%
  group_by(gene_id) %>%
  summarise(transcript_count = n_distinct(transcript_id)) # n_distinct() counts the number of distinct combinations in a set of vectors.

# Store the genes with more than or equal to 2 transcripts using filter
transcripts = transcript_count_per_gene %>%
  filter(transcript_count >= 2)

# Filter the merged data frame based on the gene id and sort it, so now filtered data only has genes which have >= 2 isoforms
filtered_data = count_mapped_data %>%
  semi_join(transcripts, by = "gene_id") %>% arrange(gene_id) # semi_join() returns all rows from x with a match in y and arrange() orders the rows of a data frame by the values of selected columns.

filtered_data = merge(transcripts, filtered_data, by = "gene_id")
#72682 obsv

filtered_data = na.omit(filtered_data) # Omit any NAs in the data
#72682 obsv

save(filtered_data, file = "new_filtered_data.Rdata")

###### 
# Calculate and sort the variances across the transcripts
load("C:/Users/pooja/new_filtered_data.Rdata")
filtered_data_new = filtered_data %>% # Filter out samples that have no gene_symbol
  filter(gene_symbol != "")
#72531
ordered_data = as.data.frame(arrange(filtered_data_new, desc(transcript_count)))

#dat = ordered_data[ordered_data$gene_id == "ENSG00000091831", ]

transcript_variances = list()
transcript_median = list()
###### VARIANCE FOR ALL ########
#grep("ESR1", ordered_data$gene_symbol)
for (transcript_id in unique(ordered_data$transcript_id)) {
  # Subset data for each transcript
  dat = ordered_data[ordered_data$gene_id == transcript_id, ]
  
  # Extract count matrix excluding non-numeric columns
  transcript_data_count_matrix = as.matrix(dat[, -c(1:9)])
  
  # Calculate row variances
  row_variances = apply(transcript_data_count_matrix, 1, sd)
  row_median = apply(transcript_data_count_matrix,1, median)
  # Append variances to transcript_variances list
  transcript_variances = append(transcript_variances, list(row_variances))
  transcript_median = append(transcript_median, list(row_median))
}
ordered_data$transcript_variance = unlist(transcript_variances)
ordered_data$transcript_median = unlist(transcript_median)
############ X #################

# Re-arrange the dataframe
variance_data = ordered_data %>%
  select(gene_id, transcript_count, gene_symbol, transcript_id, transcript_variance, transcript_median, transcript_start, transcript_end, transcript_start_site, everything())

#save(variance_data, file = "variance_data.Rdata")
###### TEST SAMPLE #######
od = ordered_data[184:500,]
for (transcript_id in unique(od$transcript_id)) {
  # Subset data for each transcript
  dat = od[od$transcript_id == transcript_id, ]

  # Extract count matrix excluding non-numeric columns
  transcript_data_count_matrix = as.matrix(dat[, -c(1:9)])
  
  # Calculate row variances
  row_variances = apply(transcript_data_count_matrix, 1, sd)
  row_median = apply(transcript_data_count_matrix,1, median)
  # Append variances to transcript_variances list
  transcript_variances = append(transcript_variances, list(row_variances))
  transcript_median = append(transcript_median, list(row_median))
}

# Append transcript_variances back to data
od$transcript_variance = unlist(transcript_variances)
od$transcript_median = unlist(transcript_median)
############ X #################

# Re-arrange the dataframe
variance_data = od %>%
  select(gene_id, transcript_count, gene_symbol, transcript_id, transcript_variance, transcript_median, transcript_start, transcript_end, transcript_start_site, everything())

#save(variance_data, file = "variance_data.Rdata")
#load("D:/BINP37/variance_data.Rdata")

###### GENE RANK ########
sorted_data_test = as.data.frame(arrange(variance_data, desc(transcript_variance)))

###### SORT PER GENE #######
# Sort the data, so that each gene is arranged in descending order of the transcript variances and remove transcripts with sd = 0
sorted_data = list()
for(i in unique(variance_data$gene_id)) {
  gene_info = variance_data[variance_data$gene_id == i,]
  arranged_info = gene_info %>%
    arrange(desc(transcript_variance)) 
  sorted_data = rbind(sorted_data, arranged_info)
}

combined_list = list()

# Determine the number of intervals
num_intervals = ceiling(length(sorted_data) / ncol(variance_data))

# Loop through each interval
for (i in 1:num_intervals) {
  # Determine the start and end indices for the current interval
  start_index = (i - 1) * ncol(variance_data) + 1
  end_index = min(i * ncol(variance_data), length(sorted_data))
  
  # Extract data for the current interval
  interval_data = sorted_data[start_index:end_index]
  
  # Store the data in the combined list
  combined_list[[i]] = interval_data
}

# Combine all intervals column-wise
variance_arranged_data = as.data.frame(do.call(cbind, combined_list))

#load("D:/BINP37/variance_arranged_data.Rdata")
# X #

###### FILTER BASED ON SORTED PER GENE ##########
# filter out sum zeros
variance_arranged_data_filtered = variance_arranged_data %>%
  group_by(gene_id) %>%
  filter(transcript_variance > 1)
# Filter out rows with 5% or 10% samples which have less than 1 FPKM
for(i in unique(variance_arranged_data_filtered$gene_id)){
  dat = variance_arranged_data_filtered[variance_arranged_data_filtered$gene_id == i,]
  dat_matrix = as.matrix(dat[, -c(1:11)])
}

rownames(dat_matrix) = rep("FPKM", nrow(dat_matrix))
num_samples = ncol(dat_matrix)

# Count the number of samples with FPKM more than 1 for each row
above_1_FPKM = rowSums(dat_matrix > 1, na.rm = TRUE)

# Calculate the percentage of samples with FPKM more than 1 for each row
percent_val = above_1_FPKM / num_samples

# Filter rows where less than 5% of the samples have FPKM less than 1
filtered_dat_matrix = dat_matrix[percent_val >= 0.05, , drop = FALSE]

#save(variance_arranged_data_filtered, file = "variance_arranged_data_filtered.Rdata")
#load("D:/BINP37/variance_arranged_data_filtered.Rdata")


########## FILTER BASED ON GENE RANK #############
variance_filtered = sorted_data_test %>%
  group_by(gene_id) %>%
  filter(transcript_variance > 1)
# Calculate the percentage of samples that have > 1 FPKM
percent_list = list()
for(i in unique(variance_filtered$transcript_id)){
  dat = variance_filtered[variance_filtered$transcript_id == i,]
  dat_matrix = as.matrix(dat[, -c(1:11)])
  rownames(dat_matrix) = c(dat$transcript_id)
  num_samples = ncol(dat_matrix)
  above_1_FPKM = rowSums(dat_matrix > 1, na.rm = TRUE)
  percent_val = (above_1_FPKM / num_samples) * 100
  percent_list = append(percent_list, percent_val)
}

variance_filtered$percent_samples = unlist(percent_list)
vd = variance_filtered %>%
  select(gene_id, transcript_count, gene_symbol, transcript_id, transcript_variance, transcript_median, transcript_start, transcript_end, transcript_start_site, percent_samples, everything())
# X #



######
library(ComplexHeatmap)

#hc = hclust(dist((gene_matrix)), method = "ward.D2")
#column_order = colnames(gene_matrix)[hc$order]

gene_data = filtered_data[filtered_data$gene_symbol =="ESR1",]
# Remove the first two columns
gene_count_only = gene_data[,-c(1:9)]

# Convert the remaining data frame to a matrix
gene_matrix = as.matrix(gene_count_only)

# Rename columns and rows
colnames(gene_matrix) = sapply(sample_name_split(colnames(gene_matrix)), `[`, 1)
rownames(gene_matrix) = c(gene_data$transcript_id)

# Define your colour pallete
cols_exp = colorRampPalette(c("red", 'white', "green"))(n=20)

# Colour list 
col_list_NCN = c("Basal"= "red", "Her2"= "purple", "LumA"= "blue", "LumB" = "lightblue", "Normal"= "green","unclassified"= "gray")
col_list_NHG = c("1" = "gray",  "2" = "brown" ,"3" = "darkred")
col_list_DRFI = c("0" = "pink","1" = "violet")
# Create a categorical variable, and order of levels
NCN_levels = factor(isoform_info$PAM50_NCN, levels = names(col_list_NCN)) 

NHG_levels = factor(isoform_info$NHG, levels = names(col_list_NHG))

DRFI_levels = factor(isoform_info$DRFI_bin, levels = names(col_list_DRFI))

# Create  a vector which holds the associated colors for each catergory
#col_vec = col_list_NCN[isoform_info$PAM50_NCN[match(colnames(gene_matrix), isoform_info$GEX_assay)]]

# Define a column annotation for Heatmap
column_annotation = HeatmapAnnotation(
  NCN = NCN_levels,
  NHG = NHG_levels, 
  DRFI = DRFI_levels, 
  col = list(NCN = col_list_NCN, NHG = col_list_NHG, DRFI = col_list_DRFI), 
  na_col = "beige"
  )

#pdf(paste("D:/BINP37/Project/Figures/", "KRIT1_heatmap_expression_test1.pdf", sep = ""))
h1 = (Heatmap(gene_matrix,
        col = cols_exp, 
        name ="Expression levels", 
        column_title = "Expression of ESR1 gene",
        cluster_columns = TRUE,
        show_heatmap_legend = TRUE,
        clustering_method_columns = "ward.D2",
        column_title_gp = gpar(fontsize = 10),
        column_names_gp =  gpar(fontsize = 0.5),
        row_names_gp = gpar(fontsize = 5),
        top_annotation = column_annotation,
        width = 12))
  dev.off()


# Draw the second heatmap  
#pdf(paste("D:/BINP37/Project/", "KRIT1_heatmap_methylation_test1.pdf", sep = ""))


h2 = (Heatmap(beta.data,
             col = my_palette, 
             name ="Methylation degree",
             cluster_rows = TRUE,
             column_order = column_order(h1),
             clustering_method_columns = "ward.D2",
             column_title = paste(my.gene,"promoter: upS=",my.upS,"dnS=",my.dnS),
             show_heatmap_legend = TRUE,
             column_title_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 0.5),
             row_names_gp = gpar(fontsize = 5),
             width = 12))

dev.off()

identical(column_order(h1),column_order(h2))

# Plot both heatmaps together
library(gridExtra)
#library(plotly)
# Arrange the heatmaps in a grid
grob_h1 = ggplotify::as.grob(h1)
grob_h2 = ggplotify::as.grob(h2)

# Arrange the grobs in a grid
grid.arrange(grob_h1, grob_h2, ncol = 1)

# library(ComplexHeatmap)

# First heatmap
h1 <- Heatmap(gene_matrix,
              col = my_pallete, 
              name = "Expression levels", 
              column_title = "Expression of ESR1 gene",
              cluster_columns = TRUE,
              show_heatmap_legend = TRUE,
              clustering_method_columns = "ward.D2",
              column_title_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 0.5),
              row_names_gp = gpar(fontsize = 5),
              top_annotation = column_annotation,
              width = 12)

# Second heatmap
h2 <- Heatmap(beta.data,
              col = my_palette, 
              name = "Methylation degree",
              cluster_rows = FALSE,
              clustering_method_columns = "ward.D2",
              column_order = column_order(h1),
              column_title = paste(my.gene, "promoter: upS=", my.upS, "dnS=", my.dnS),
              show_heatmap_legend = TRUE,
              column_title_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 0.5),
              row_names_gp = gpar(fontsize = 5),
              width = 12)


##### 
# box plot

library(ggplot2)
library(reshape2)
library(gridExtra)


gene_data = filtered_data[filtered_data$gene_symbol == "ESR1", ]
# Remove the first two columns
gene_count_only = gene_data[,-c(1:9)]

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

# Extract unique ENST_IDs
unique_ENST = unique(melted_data$ENST_ID)

# List to store plots
plots_list = list()
#colors_list <- list()
# Loop through each ENST_ID and create a plot
for (enst in unique_ENST) {
  # Subset data for the current ENST_ID
  enst_data = subset(melted_data, ENST_ID == enst)
  
  # Create plot for the current ENST_ID
  p = ggplot(enst_data, aes(x = factor(NCN), y = Count, fill = factor(NCN))) +
    geom_boxplot() +
    labs(x = "NCN Levels", y = "Count", title = paste("ENST ID:", enst, "  ")) +
    theme_minimal()
  
  # Store the plot in the list
  plots_list[[enst]] <- p
  #colors_list[[enst]] <- ggplot_build(p)$data[[1]]$fill
}

# Combine all plots into one
combined_plot = do.call(grid.arrange, c(plots_list, ncol = 3))  # Adjust ncol as needed

ggsave("combined_plot.pdf", combined_plot, width = 20, height = 15)

#for (enst in unique_ENST) {
#  cat("ENST ID:", enst, "Colors:", unique(colors_list[[enst]]), "\n")
#}


######

install.packages("biomartr")
library(biomartr)
getGTF(
  db = "ensembl",
  organism = "Homo sapiens",
  remove_annotation_outliers = FALSE,
  path = file.path("ensembl", "annotation"),
  release = NULL,
  mute_citation = FALSE
)



grep("S006452", names(column_order(h1)))
colnames(gene_matrix)[476]
