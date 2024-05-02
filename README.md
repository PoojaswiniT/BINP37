# Project description
Visualising the correlation of Isoform expression with DNA methylation in breast cancer cells. 

## Create a mapping of ENSTxxx to ENSG with gene symbol for each isoform
### Files used
- transcription_count_matrix.csv

### Methods
- In R studio, load the transcript_count_matrix.csv file.
- The ENST ids have versions, so split the IDs by "." to match the BioMart naming.
- Using biomart package, obtain all the atributes- "ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "transcript_start", "transcript_end", "transcript_start_site".
- Then, using merge function with dplyr, merge the transcription_count_matrix with transcript_to_gene_mapping dataframe based on ensembl_transcript_id.

```R
# Read the data
tm = read.csv("transcript_count.csv", sep = ",")

#######
#Splitting the ID column 

split_ID = function(transcript_matrix){
  # Split the ID based on dot
  id_split = strsplit(as.character(transcript_matrix$transcript_id), "\\.")
  
  # Extract values before and after the dot
  transcript_matrix$transcript_id = sapply(id_split, `[`, 1)
  transcript_matrix$version = sapply(id_split, `[`, 2)
  
  # Store the read counts as well
  transcript_matrix = transcript_matrix %>%
    dplyr::select(transcript_id, version, everything())
  return(transcript_matrix)
}

transcript_count_matrix = split_ID(tm)

#####
#Mapping ENSTxxx to ENSGxxx, gene symbol and other attributes
# Load the required libraries
library(biomaRt)
library(dplyr)

# Connect to the Ensembl database
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# List available attributes
#listAttributes(ensembl)

# Define the attributes to retrieve
attributes1 = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name")
attributes2 = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name")
#attributes3 = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name")

transcript_to_gene_mapping1 = getBM(attributes = attributes1, mart = ensembl)
transcript_to_gene_mapping2 = getBM(attributes = attributes2, mart = ensembl)
#transcript_to_gene_mapping2 = getBM(attributes = attributes3, mart = ensembl)

enst_to_ensg_mapping = merge(transcript_to_gene_mapping1, transcript_to_gene_mapping2, by = "ensembl_transcript_id")

# Change the column name in accordance to the ensembl format
colnames(transcript_count_matrix)[1] = "ensembl_transcript_id"

# Merge the data frames based on the ensembl_transcript_id column
mapped_data = merge(transcript_count_matrix, enst_to_ensg_mapping, by = "ensembl_transcript_id", all.x = TRUE)
# all.x = TRUE : all rows from the left dataframe are kept and matching rows from the right dataframe are included, with NA values filled in for non-matching rows

# Rearrange the columns
mapped_data = merged_data %>%
  select(ensembl_transcript_id, version, ensembl_gene_id, external_gene_name, everything())

#Matching the annotations from the EPIC data 

load("Patient_annotations/Annotations_matched_EPIC_ERpHER2n.RData") # Load the matched patients data

# Obtain all the annotations
isoform_info = data.frame(
  GEX_assay = annotations$GEX.assay, # Sample names
  PAM50_NCN = annotations$PAM50_NCN, # Cancer type
  NHG = annotations$NHG, # Nottingham Histological Grade (NHG)
  DRFI_bin = annotations$DRFIbin # Distant-recurrence free interval (DRFI) probabilities
  )

# Merge data frames based on the sample columns
columns = list()

for (i in isoform_info$GEX_assay){ # loop through all the sample names
  # If the sample name in the 515 set matches the sample name in the count matrix
  if (i %in% setdiff(colnames(mapped_data), c("ensembl_gene_id", "ensembl_transcript_id", "version", "external_gene_name", "chromosome_name", "strand", "transcript_start", "transcript_end", "transcript_start_site"))){
    columns[[i]] = mapped_data[, i] # Store the sample columns that match     
  }
}

new_mapped_data = as.data.frame(do.call(cbind,columns)) # column bind the list

# add the remaining info
new_mapped_data$gene_id = mapped_data$ensembl_gene_id
new_mapped_data$transcript_id = mapped_data$ensembl_transcript_id
new_mapped_data$gene_symbol = mapped_data$external_gene_name
new_mapped_data$chromosome = mapped_data$chromosome_name
new_mapped_data$strand = mapped_data$strand
new_mapped_data$transcript_start = mapped_data$transcript_start
new_mapped_data$transcript_end = mapped_data$transcript_end
new_mapped_data$transcript_start_site = mapped_data$transcription_start_site

# Rearrange the data frame
count_mapped_data = new_mapped_data %>%
  dplyr::select(gene_id, gene_symbol, transcript_id, chromosome, strand, transcript_start, transcript_end, transcript_start_site, everything())


#####
# Counting number of transcripts for each gene
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
#72682 obsv #18232

filtered_data = na.omit(filtered_data) # Omit any NAs in the data
#72682 obsv #18232

#####
#Gene_expression with Complex Heatmap
load("D:/BINP37/filtered_data.Rdata")

# Arrange the filtered data by descending order of the number of transcripts
ordered_data = as.data.frame(arrange(filtered_data, desc(transcript_count)))

gene_data = data.frame()

isoform_info$GEX_assay = sapply(sample_name_split(annotations$GEX.assay), `[`, 1) # Simplifying the sample names
col_list = c("Basal"= "red", "Her2"= "purple", "LumA"= "blue", "LumB" = "lightblue", "Normal"= "green","unclassified"= "gray") # Standard colours 
# Gene expression with complexc heatmap
for(j in unique(ordered_data$gene_symbol)) {
  
  gene_data = ordered_data[ordered_data$gene_symbol == j, ]
  gene_data = gene_data[,-c(1:2,5:9)]
  
  # Remove the first two columns
  gene_count_only = gene_data[,-c(1:2)]
  
  # Convert the remaining data frame to a matrix
  gene_matrix = as.matrix(gene_count_only)
  
  # Rename columns and rows
  #colnames(gene_matrix) = sapply(sample_name_split(colnames(gene_matrix)), `[`, 1)
  #rownames(gene_matrix) = c(gene_data$transcript_id)
  
  
  # Define the colour pallete
  my_pallete = colorRampPalette(c("red", 'white', "green"))(n=20)
  
  # Colour list 
  col_list_NCN = c("Basal"= "red", "Her2"= "purple", "LumA"= "blue", "LumB" = "lightblue", "Normal"= "green","unclassified"= "gray")
  col_list_NHG = c("1" = "gray",  "2" = "brown" ,"3" = "darkred")
  col_list_DRFI = c("0" = "pink","1" = "violet")
  # Create a categorical variable, and order of levels
  NCN_levels = factor(isoform_info$PAM50_NCN, levels = names(col_list_NCN)) 
  
  NHG_levels = factor(isoform_info$NHG, levels = names(col_list_NHG))
  
  DRFI_levels = factor(isoform_info$DRFI_bin, levels = names(col_list_DRFI))
  
  # Define a column annotation for Heatmap
  column_annotation = HeatmapAnnotation(
    NCN = NCN_levels,
    NHG = NHG_levels, 
    DRFI = DRFI_levels, 
    col = list(NCN = col_list_NCN, NHG = col_list_NHG, DRFI = col_list_DRFI), 
    na_col = "beige"
  )
  
  pdf(paste("D:/BINP37/Project/Complex_heatmaps/", j, "_heatmap.pdf", sep = ""))
  draw(Heatmap(gene_matrix,
          col = my_pallete, 
          name ="Expression levels", 
          column_title = paste0("Expression of ",j," gene"), 
          show_heatmap_legend = TRUE, 
          #column_order = colnames(gene_matrix),
          column_title_gp = gpar(fontsize = 10),
          column_names_gp =  gpar(fontsize = 0.5),
          row_names_gp = gpar(fontsize = 5),
          top_annotation = column_annotation,
          width = 12))
  #dev.off()
}