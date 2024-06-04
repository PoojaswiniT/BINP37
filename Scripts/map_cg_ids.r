##### Mapping the respective CpG ids from beta matrix ####

#### Load the library ####
library(GenomicRanges)

#Load the annotations to work with
load("D:/BINP37/Pooja/tracks/Patient_annotations/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData")


#### Subselect FU=TRUE & ERpHER2n ##
ii = which( (pam50.frame$ClinGroup=="ERpHER2nLNn") | (pam50.frame$ClinGroup=="ERpHER2nLNp") )
pam50.frame=pam50.frame[ii,]

ii=which( pam50.frame$Follow.up.cohort  )
pam50.frame=pam50.frame[ii,]
#4924 samples
###

# Load the promoter annotation object for EPIC data
load("D:/BINP37/Pooja/getPromoter.RData")

# Load the function to get probes
load("D:/BINP37/New data/get_cg.Rdata")

### Define the genes to analyze ####
load("variance_percent.Rdata")
my.genes= c(vd$gene_symbol)

### Load the BASE-processed EPIC data, edit names, match to WGS ###
load("D:/BINP37/Given data/scanb_base.RData")
nn=colnames(scanb_base)
nn=gsub(".l.d.mth","",nn,fixed=T)
nn=gsub(".l2.d.mth","",nn,fixed=T)
colnames(scanb_base)=nn
rownames(scanb_base)=scanb_base$illuminaID
beta.matrix=as.matrix(scanb_base[,-1])
rownames(beta.matrix)=scanb_base$illuminaID
rm(scanb_base)
gc();
ii=match(pam50.frame$Sample,colnames(beta.matrix))
ii.i=which(!is.na(ii))
beta.matrix=beta.matrix[,ii[ii.i]]
annotations=pam50.frame[ii.i,]
stopifnot(identical(annotations$Sample,colnames(beta.matrix)))
gc()
print(dim(beta.matrix))
###

### Iterate over genes ###
my.upS=10000
my.dnS=2000

gene_info_list = list()

# Loop through each unique gene
for (i in unique(vd$gene_symbol)) {
    gene_data = vd[vd$gene_symbol == i,]
    
    if (nrow(gene_data) == 0) {
        print(paste("no match found for:", i))
        next
    }
    
    if (gene_data$strand[1] == "1") {
        strand = "+"
        gene_data$strand = "+"
    } else {
        strand = "-"
        gene_data$strand = "-"
    }
    
    my.upS = 10000
    my.dnS = 2000
    myProm = as.data.frame(get_cg(gene_data, my.upS = my.upS, my.dnS = my.dnS))
    myProm = myProm[!duplicated(myProm), ]
    
    if (nrow(myProm) == 0) {
        print(paste("no match", i, ":", i))
        next
    }
    
    # Collect the actual beta data
    ii = match(myProm$probes, rownames(beta.matrix))
    
    # Check if 'myProm$probes' has valid matches in 'beta.matrix'
    if (all(is.na(ii))) {
        print(paste("no valid beta data for:", i))
        next
    }
    
    # Collect gene information
    gene_info = list(
        gene_name = i,
        ensg = gene_data$gene_id[1],
        cpgs = myProm$probes,
        strand = gene_data$strand[1]
    )
    
    # Append the gene information to the list
    gene_info_list[[length(gene_info_list) + 1]] = gene_info
}

# Convert the list to a data frame
gene_info_df = do.call(rbind, lapply(gene_info_list, as.data.frame))

# View the data frame
View(gene_info_df)

# Convert list columns to appropriate data types
cg_ids = data.frame(
    gene_name = unlist(gene_info_df$gene_name),
    ensg = unlist(gene_info_df$ensg),
    cpgs = unlist(gene_info_df$cpgs),
    strand = unlist(gene_info_df$strand))
