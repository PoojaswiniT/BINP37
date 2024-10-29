###

# Description - the following code reads the given transcript count matrix and maps the transcript
# ids to their respective gene ids, and other attritributes using biomaRt. 
# Then it is matched with annotations of the EPIC data provided.


##### Load libraries
library(dplyr)
library(tidyr)
library(biomaRt)

#####
# Read the transcript count matrix data
tm = read.csv("/Users/pooja/Documents/BINP37/transcript_fpkm_matrix-4.3.csv", sep = "\t")
colnames(tm)[1] = "transcript_id" #rename the column

#Splitting the version number of transcript IDs to match the biomart usage
source("/Users/pooja/Documents/BINP37/split_function.R")
transcript_count_matrix = as.data.frame(split_ID(tm))

#### Mapping ENSG IDs and other attributes #####

#Connect to the Ensembl database
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#List available attributes
list_attributes  = as.data.frame(listAttributes(ensembl))

#Store the attributes, separately as too many variables at once leads to timeout error

attributes1 = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "chromosome_name", "strand")
attributes2 = c("ensembl_transcript_id", "transcript_start", "transcript_end", "transcription_start_site")

#Get the attributes from the package
transcript_to_gene_mapping1 = getBM(attributes = attributes1, mart = ensembl)
transcript_to_gene_mapping2 = getBM(attributes = attributes2, mart = ensembl)

#Merge all the information 
enst_to_ensg_map = merge(transcript_to_gene_mapping1, transcript_to_gene_mapping2, by = "ensembl_transcript_id", all = TRUE)

#Change the column name of the input matrix in accordance to the ensembl format
colnames(transcript_count_matrix)[1] = "ensembl_transcript_id"

#Merge the data frames based on the ensembl_transcript_id column
mapped_data = merge(transcript_count_matrix, enst_to_ensg_map, by = "ensembl_transcript_id", all.x = TRUE)

#Rearrange the columns
mapped_data = mapped_data %>%
  dplyr::select(ensembl_gene_id, ensembl_transcript_id, version, external_gene_name, chromosome_name, strand, transcript_start, transcript_end, transcription_start_site, everything())
#80967

#Omit all the NAs
mapped_data = na.omit(mapped_data)
#77526

#save(mapped_data, file = "Documents/BINP37/new_mapped_data.Rdata")


##### Matching the annotations from the EPIC data #####
#load("D:/BINP37/new_mapped_data.Rdata") # Mapped data

load("/Users/pooja/Documents/BINP37/Annotations_matched_EPIC_ERpHER2n.RData") # the annotation tracks provided

md = mapped_data # For working

isoform_info = data.frame(
  GEX_assay = annotations$GEX.assay, # Sample names
  PAM50_NCN = annotations$PAM50_NCN, # Cancer type
  NHG = annotations$NHG, # Nottingham Histological Grade (NHG)
  DRFI_bin = annotations$DRFIbin # Distant-recurrence free interval (DRFI) probabilities
)

# Merge data frames based on the sample columns
columns = list() # Initialize a list

for (i in isoform_info$GEX_assay){ # Loop through all the sample names
  # If the sample name in the 515 set matches the sample name in the count matrix
  if (i %in% setdiff(colnames(md), c("ensembl_gene_id", "ensembl_transcript_id", "version", "external_gene_name", "chromosome_name", "strand", "transcript_start", "transcript_end", "transcript_start_site"))){
    columns[[i]] = md[, i] # Store the sample columns that match
  }
}

new_mapped_data = as.data.frame(do.call(cbind,columns)) # Column bind all the values to a dataframe

#Add the remaining information
new_mapped_data$gene_id = mapped_data$ensembl_gene_id
new_mapped_data$transcript_id = mapped_data$ensembl_transcript_id
new_mapped_data$gene_symbol = mapped_data$external_gene_name
new_mapped_data$chromosome = mapped_data$chromosome_name
new_mapped_data$strand = mapped_data$strand
new_mapped_data$transcript_start = mapped_data$transcript_start
new_mapped_data$transcript_end = mapped_data$transcript_end
new_mapped_data$transcript_start_site = mapped_data$transcription_start_site

#Re-arrange the data frame
count_mapped_data = new_mapped_data %>%
  dplyr::select(gene_id, gene_symbol, transcript_id, chromosome, strand, transcript_start, transcript_end, transcript_start_site, everything())
#77526 rows # 523 columns

#save(count_mapped_data, file = "count_mapped_data.Rdata")
