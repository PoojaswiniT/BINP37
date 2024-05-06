#####
# Load libraries
library(dplyr)
library(tidyr)

#####
# Read the transcript count matrix data
tm = read.csv("D:/BINP37/New data/transcript_fpkm_matrix-4.3.csv", sep = "\t")

colnames(tm)[1] = "transcript_id"
#######
#Splitting the versions of transcript IDs
load("D:/BINP37/New data/split_function.R")
transcript_count_matrix = as.data.frame(split_ID(tm))
#####
# Mapping ENSG IDs and other attributes
library(biomaRt)

# Connect to the Ensembl database
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Package version of biomart
#packageVersion("biomaRt")

# List of datasets available
#listDatasets(mart = ensembl)

# Available attributes
list_attributes  = as.data.frame(listAttributes(ensembl))

# Store the attributes required separately as too many variables at once leads to timeout error

attributes1 = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "chromosome_name", "strand")
attributes2 = c("ensembl_transcript_id", "transcript_start", "transcript_end", "transcription_start_site")

# Get the attributes from the package
transcript_to_gene_mapping1 = getBM(attributes = attributes1, mart = ensembl)
transcript_to_gene_mapping2 = getBM(attributes = attributes2, mart = ensembl)

# Merge all the information 
enst_to_ensg_map = merge(transcript_to_gene_mapping1, transcript_to_gene_mapping2, by = "ensembl_transcript_id", all = TRUE)

# Change the column name of the input matrix in accordance to the ensembl format
colnames(transcript_count_matrix)[1] = "ensembl_transcript_id"

# Merge the data frames based on the ensembl_transcript_id column
mapped_data = merge(transcript_count_matrix, enst_to_ensg_map, by = "ensembl_transcript_id", all.x = TRUE)

# Rearrange the columns
mapped_data = mapped_data %>%
  dplyr::select(ensembl_gene_id, ensembl_transcript_id, version, external_gene_name, chromosome_name, strand, transcript_start, transcript_end, transcription_start_site, everything())
#80967

# Omit all the NAs
mapped_data = na.omit(mapped_data)
# 77526
#save(mapped_data, file = "new_mapped_data.Rdata")

# Use dplyr::select when select doesn't work

#####
# Matching the annotations from the EPIC data 
load("D:/BINP37/new_mapped_data.Rdata") # Mapped data

load("D:/BINP37/Pooja/tracks/Patient_annotations/Annotations_matched_EPIC_ERpHER2n.RData") # the annotation tracks provided

load("D:/BINP37/sample_name_split.R") # A function to split the sample names


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

# Add the remaining information
new_mapped_data$gene_id = mapped_data$ensembl_gene_id
new_mapped_data$transcript_id = mapped_data$ensembl_transcript_id
new_mapped_data$gene_symbol = mapped_data$external_gene_name
new_mapped_data$chromosome = mapped_data$chromosome_name
new_mapped_data$strand = mapped_data$strand
new_mapped_data$transcript_start = mapped_data$transcript_start
new_mapped_data$transcript_end = mapped_data$transcript_end
new_mapped_data$transcript_start_site = mapped_data$transcription_start_site

# Re-arrange the data frame
count_mapped_data = new_mapped_data %>%
  dplyr::select(gene_id, gene_symbol, transcript_id, chromosome, strand, transcript_start, transcript_end, transcript_start_site, everything())
# 77526 rows # 523 columns


# Mapping the respective CpG ids
library(GenomicRanges)

# Load the annotations to work with
load("D:/BINP37/Pooja/tracks/Patient_annotations/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData")
############

#### Subselect FU=TRUE & ERpHER2n ##
ii <- which( (pam50.frame$ClinGroup=="ERpHER2nLNn") | (pam50.frame$ClinGroup=="ERpHER2nLNp") )
pam50.frame<-pam50.frame[ii,]

ii<-which( pam50.frame$Follow.up.cohort  )
pam50.frame<-pam50.frame[ii,]
#4924 samples
###

### Load the promoter annotation object for EPIC data ###
load("D:/BINP37/Pooja/getPromoter.RData")
###
### Define the genes to analyze ####
my.genes<- c(count_mapped_data$gene_symbol)
###

assign('.Last',  function() {system('R')}, envir = globalenv())
quit(save = 'yes')

### Load the BASE-processed EPIC data, edit names, match to WGS ###
load("D:/BINP37/Given data/scanb_base.RData")
nn<-colnames(scanb_base)
nn<-gsub(".l.d.mth","",nn,fixed=T)
nn<-gsub(".l2.d.mth","",nn,fixed=T)
colnames(scanb_base)<-nn
rownames(scanb_base)<-scanb_base$illuminaID
beta.matrix<-as.matrix(scanb_base[,-1])
rownames(beta.matrix)<-scanb_base$illuminaID
rm(scanb_base)
gc();
ii<-match(pam50.frame$Sample,colnames(beta.matrix))
ii.i<-which(!is.na(ii))
beta.matrix<-beta.matrix[,ii[ii.i]]
annotations<-pam50.frame[ii.i,]
stopifnot(identical(annotations$Sample,colnames(beta.matrix)))
gc()
print(dim(beta.matrix))
###


##Get promoter data ###
ls()
#[1] "geneCoords" "getProm"    "probeAnno" 


### Iterate over genes - save to a list with each entry a gene in the form of a list. names of list = gene name ###
gene.set.list<-list()
nn.gene<-c()#to hold names
my.upS<-10000
my.dnS<-2000


for(i in 1:length(my.genes[1:10])){
  my.gene<-my.genes[i]
  #print(paste(i,":",my.gene))
  myProm<-getProm(symbol=my.gene,upS=my.upS,dnS=my.dnS)
  if(length(myProm)==0){
    print(paste("no match",i,":",my.gene))
    next
    #no match
  }else{
    #collect CpG info data
    
    #collect the actual beta data
    ii<-match(myProm$probes,rownames(beta.matrix))
    ii.i<-which(!is.na(ii))
    #beta.data<-beta.matrix[ii[ii.i],]
    
    #safety check
    #if(nrow(beta.data)==0) next
    #
    
    nn.gene<-c(nn.gene,my.gene)
    #myProm$beta<-beta.data
    gene.set.list[[length(gene.set.list)+1]]<-myProm
  }
}
names(gene.set.list)<-nn.gene
cg_data =  bind_rows(gene.set.list)

###

### Save list ###
save(gene.set.list,file=paste("GenePromoterSets_upS",my.upS,"_dnS",my.dnS,".RData",sep=""))
load("D:/BINP37/New data/GenePromoterSets_upS10000_dnS2000.RData")
####

# Split the gene based on the pipe
gene_split = function(cg){
  g_split = strsplit(as.character(cg$gene), "\\|")
  cg$gene = sapply(g_split, `[`, 2) # Since the gene symbol is on the right 
  cg = cg%>%
    dplyr::select(gene, everything())
  return(cg)
  }

cg_new = as.data.frame((gene_split(cg_data)))
colnames(cg_new)[1] = "gene_symbol"
colnames(count_mapped_data)[2]
cg_merged = merge(count_mapped_data, cg_new, by= "gene_symbol")


cg_data_arranged = cg_merged %>%
dplyr::select(gene_id, gene_symbol, transcript_id, strand, chromosome, chrom, transcript_start, transcript_end, transcript_start_site, tssStart, tssDist, probes, probeStart, everything())

save(cg_data_arranged, file = "cg_data.Rdata")
load("D:/BINP37/cg_data.Rdata")
