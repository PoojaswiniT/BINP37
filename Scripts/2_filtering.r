###

# Description - the code filters the genes which have less than 2 transcripts, as they might not provide much information;
# also arranges the genes in a ranking order based on the transcript variance across rows.


#### Load libraries ####
library(dplyr)
library(tidyr)

#Count number of transcripts
#Group by gene id, count the number of distinct transcript ids for each gene
transcript_count_per_gene = count_mapped_data %>%
  group_by(gene_id) %>%
    summarise(transcript_count = n_distinct(transcript_id)) # n_distinct() counts the number of distinct combinations in a set of vectors.

#Store the genes with more than or equal to 2 transcripts
transcripts = transcript_count_per_gene %>%
  filter(transcript_count >= 2)

#Filter the merged data frame based on the gene id and sort it, so now filtered data only has genes which have >= 2 isoforms
filtered_data = count_mapped_data %>%
  semi_join(transcripts, by = "gene_id") %>% arrange(gene_id) 

#semi_join() returns all rows from x with a match in y and arrange() orders the rows of a data frame by the values of selected columns.

filtered_data = merge(transcripts, filtered_data, by = "gene_id")
#72682 obsv

filtered_data = na.omit(filtered_data) # Omit any NAs in the data
#72682 obsv

#save(filtered_data, file = "new_filtered_data.Rdata")

#### Calculate and sort the variances across the transcripts ####
#load("C:/Users/pooja/new_filtered_data.Rdata")

filtered_data_new = filtered_data %>% # Filter out samples that have no gene_symbol
  filter(gene_symbol != "")
#72531

ordered_data = as.data.frame(arrange(filtered_data_new, desc(transcript_count))) # not really required but good to visualise

transcript_variances = list()
transcript_median = list()

for (transcript_id in unique(ordered_data$transcript_id)) {
  #Subset data for each transcript
  dat = ordered_data[ordered_data$transcript_id == transcript_id, ]
  
  #Extract count matrix with only numeric values
  transcript_data_count_matrix = as.matrix(dat[, -c(1:9)])
  
  #Calculate row variances using standard deviation
  row_variances = apply(transcript_data_count_matrix, 1, sd)
  row_median = apply(transcript_data_count_matrix,1, median) #Calculate median expression
  #Append variances to transcript_variances list
  transcript_variances = append(transcript_variances, list(row_variances))
  transcript_median = append(transcript_median, list(row_median))
}

ordered_data$transcript_variance = unlist(transcript_variances)
ordered_data$transcript_median = unlist(transcript_median)

#Re-arrange the dataframe
variance_data = ordered_data %>%
  select(gene_id, transcript_count, gene_symbol, transcript_id, transcript_variance, transcript_median, transcript_start, transcript_end, transcript_start_site, everything())

#save(variance_data, file = "variance_data.Rdata")

###### GENE RANK ########
sorted_data_test = as.data.frame(arrange(variance_data, desc(transcript_variance))) #Sort the data frame in the descending order of transcript variance

########## FILTER BASED ON GENE RANK #############
variance_filtered = sorted_data_test %>%
  group_by(gene_id) %>%
  filter(transcript_variance > 10) #Filteer out the transcripts with value less than 10

#Calculate the percentage of samples that have > 1 FPKM
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
# Filtering samples with less than 10 percent variance
vd = variance_filtered %>%
  group_by(gene_id) %>%
  filter(percent_samples > 10) %>%
  select(gene_id, transcript_count, gene_symbol, transcript_id, transcript_variance, transcript_median, transcript_start, transcript_end, transcript_start_site, percent_samples, everything())

#save(vd, file = "variance_percent.Rdata")
#load("variance_percent.Rdata")
# variance_percent.Rdata is used in the map_cg_ids.R script