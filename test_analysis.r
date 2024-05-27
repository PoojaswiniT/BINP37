library(dplyr)
library(tidyr)
load("D:/BINP37/Pooja/tracks/Patient_annotations/Annotations_matched_EPIC_ERpHER2n.RData") # the annotation tracks provided

load("D:/BINP37/sample_name_split.R") # A function to split the sample names

isoform_info = data.frame(
  Samples = annotations$Sample, # Sample names
  PAM50_NCN = annotations$PAM50_NCN, # Cancer type
  NHG = annotations$NHG, # Nottingham Histological Grade (NHG)
  DRFI_bin = annotations$DRFIbin # Distant-recurrence free interval (DRFI) probabilities
  )


library(GenomicRanges)
library(gplots)
library(RColorBrewer)
my_palette = colorRampPalette(c("blue", "white", "yellow"))(n = 20)


load("D:/BINP37/Pooja/tracks/Patient_annotations/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData")


#### Subselect FU=TRUE & ERpHER2n ##
ii=which( (pam50.frame$ClinGroup=="ERpHER2nLNn") | (pam50.frame$ClinGroup=="ERpHER2nLNp") )
pam50.frame=pam50.frame[ii,]

ii=which( pam50.frame$Follow.up.cohort  )
pam50.frame=pam50.frame[ii,]
#4924 samples
###

### Load the promoter annotation object for EPIC data ###
load("D:/BINP37/Pooja/getPromoter.RData")


assign('.Last',  function() {system('R')}, envir = globalenv())
quit(save = 'yes')


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
#print(dim(beta.matrix))

load("D:/BINP37/New data/PAM50genes_geneAnnotation_Rel4.RData")
pam_50_genes = as.data.frame(y)
dna_repair_genes = as.data.frame(read.delim("D:/BINP37/full_pipline/DNArepairGenes_MDanderson.txt"))
load("variance_data.Rdata")


library(ComplexHeatmap)
library(grid)
library(gridExtra)
#for(i in variance_data$gene_symbol){

gene_data = variance_data[variance_data$gene_symbol == "FGFR4",]
# nrow(gene_data)
#gene_data = cg_data_arranged[cg_data_arranged$gene_symbol=="",]
if (nrow(gene_data) == 0) {
        print(paste("no match found for: ", i))
        next
}
if(gene_data$strand[1] == "1"){
    strand = "+"
    #my.upS= min(as.integer(gene_data$transcript_start_site))
    #my.dnS= max(as.integer(gene_data$transcript_start_site))
} else {
    strand = "-"
    #my.upS= max(as.integer(gene_data$transcript_start_site)) 
    #my.dnS= min(as.integer(gene_data$transcript_start_site)) 
}

#print(paste("Processing gene at index", i, "with strand", strand))
# Remove the first two columns
#gene_count_only = gene_data[,-c(1:17,533)]
gene_count_only = gene_data[, -c(1:11)]
# Convert the remaining data frame to a matrix
gene_matrix = as.matrix(gene_count_only)
# Rename columns and rows
colnames(gene_matrix) = sapply(sample_name_split(colnames(gene_matrix)), `[`, 1)
rownames(gene_matrix) = paste(gene_data$transcript_id, "\n", "TSS: ", gene_data$transcript_start_site)

cols_exp = colorRampPalette(c("red", 'white', "green"))(n=20)

gene.set.list=list()
nn.gene=c()#to hold names

my.gene = "FGFR4"

my.upS = 10000
my.dnS = 2000
load("D:/BINP37/New data/get_cg.Rdata")
myProm = get_cg(gene_data)
myProm = getProm(symbol = my.gene, upS = my.upS, dnS = my.dnS)
if(length(myProm)==0){
    print(paste("no match",i,":",my.gene))
    next
    #no match
} 


#collect the actual beta data
ii=match(myProm$probes,rownames(beta.matrix))
ii.i=which(!is.na(ii))
beta.data=beta.matrix[ii[ii.i],]

#safety check
if (!is.data.frame(beta.data) && !is.matrix(beta.data)) {
    print("beta.data is not a data frame or matrix")
    next  # Skip to the next iteration
}

if (nrow(beta.data) == 0) {
    print(paste("No beta data found for gene:", i))
    next  # Skip to the next iteration
}

#

nn.gene=c(nn.gene,my.gene)
myProm$beta=beta.data
gene.set.list[[length(gene.set.list)+1]]=myProm
bd = data.frame()
bd = beta.data
rownames(beta.data) = paste(myProm$probes," ", paste(myProm$chrom, myProm$probeStart, sep = ":"))



gene_col = rep("white", nrow(gene_matrix))
color_vector = c("black", "white")
names(color_vector) = c("black", "white")

if(nrow(gene_data) == 1){
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
            height = 10,
            column_km = 4, 
            left_annotation = HeatmapAnnotation(df = data.frame(gene_col), 
                                                        col = list(gene_col = color_vector),
                                                        which = "row", show_annotation_name = FALSE, show_legend = FALSE)))
} else{
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
            height = 10,
            column_km = 4))
}

col_order_h1 = as.integer(unlist(column_order(h1)))

tss_gene = c(unique(myProm$tssStart))
probe_start = c(unique(myProm$probeStart))
enst_tss = c(unique(gene_data$transcript_start_site))


probe_pos = numeric(length(enst_tss))
closest_probe_ids = character(length(enst_tss))
probe_col = rep("white", nrow(beta.data))

i = "FGFR4"
for (i in seq_along(enst_tss)) {
    tss_value = enst_tss[i]
    # Find the index of the closest probeStart to the current tss_gene value
    closest_probe_index = which.min(abs(probe_start - tss_value))
    # Get the closest probeStart value using the index
    closest_probe_start = probe_start[closest_probe_index]
    probe_id = unique(myProm$probes[myProm$probeStart == closest_probe_start])
     if (length(probe_id) > 1) {
        probe_id = probe_id[1]  # Take the first match if there are multiple
    }
    matching_pos = which(rownames(bd) == probe_id)
    
    if (length(matching_pos) > 0) {
        probe_pos[i] = matching_pos[1]  # Take the first match if there are multiple
        closest_probe_ids[i] = probe_id
        probe_col[probe_pos[i]] = "black"
    } else {
        warning(paste("No matching probe_id found in beta data for", probe_id))
    }
}


h2 = (Heatmap(beta.data,
            col = my_palette, 
            name = "Methylation",
            cluster_rows = FALSE,
            column_order = col_order_h1,
            clustering_method_columns = "ward.D2",
            column_title = paste(my.gene, "promoter: upS=", my.upS, "dnS=", my.dnS),
            show_heatmap_legend = TRUE,
            column_title_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 0.5),
            row_names_gp = gpar(fontsize = 5),
            width = 14,
            height = 10,
            left_annotation = HeatmapAnnotation(df = data.frame(probe_col), 
                                                col = list(probe_col = color_vector),
                                                which = "row", show_annotation_name = FALSE, show_legend = FALSE)))
col_order_h2 = column_order(h2)


# Check if column orders are identical 
stopifnot(identical((col_order_h1), (col_order_h2)))

# Arrange the heatmaps in a grid
grob_h1 = ggplotify::as.grob(h1)
grob_h2 = ggplotify::as.grob(h2)

#    pdf(paste("D:/BINP37/All_genes/",my.gene,".pdf"))
# Arrange the grobs in a grid
grid.arrange(grob_h1, grob_h2, ncol = 1)
dev.off()
#}
