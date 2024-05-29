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
load("D:/BINP37/New data/get_cg.Rdata")

library(ComplexHeatmap)
library(grid)
library(gridExtra)
for(i in dna_repair_genes$geneSymbol){
    gene_data = variance_data[variance_data$gene_symbol == i,]
    
    if (nrow(gene_data) == 0) {
            print(paste("no match found for: ", i))
            next
    }
    if(gene_data$strand[1] == "1"){
        strand = "+"
        gene_data$strand = "+"
        #my.upS= min(as.integer(gene_data$transcript_start_site))
        #my.dnS= max(as.integer(gene_data$transcript_start_site))
    } else {
        strand = "-"
        gene_data$strand = "+"
        #my.upS= max(as.integer(gene_data$transcript_start_site)) 
        #my.dnS= min(as.integer(gene_data$transcript_start_site)) 
    }
    gene_count_only = gene_data[, -c(1:11)]

    # Convert the remaining data frame to a matrix
    gene_matrix = as.matrix(gene_count_only)

    # Rename columns and rows
    colnames(gene_matrix) = sapply(sample_name_split(colnames(gene_matrix)), `[`, 1)
    rownames(gene_matrix) = paste(gene_data$transcript_id, "\n", "TSS: ", gene_data$transcript_start_site)

    cols_exp = colorRampPalette(c("red", 'white', "green"))(n=20)

    gene.set.list=list()
    nn.gene=c()#to hold names
    my.gene = i
    my.upS = 10000
    my.dnS = 2000
    myProm = as.data.frame(get_cg(gene_data, my.upS = 10000, my.dnS = 2000))

    myProm <- myProm[!duplicated(myProm), ]
    

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
    rownames(beta.data) = paste(myProm$probes," ", paste("chr", unique(gene_data$chromosome), myProm$probeStart, sep = ":"))

    gene_col = rep("#FFFFFF", nrow(gene_matrix))
    gene_color_vector = c("#000000", "#FFFFFF")
    names(gene_color_vector) = c("#000000", "#FFFFFF")

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
                width = 8, 
                height = 10,
                column_km = 4, 
                left_annotation = HeatmapAnnotation(df = data.frame(gene_col), 
                                                            col = list(gene_col = gene_color_vector),
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
                width = 8, 
                height = 10,
                column_km = 4))
    }

    col_order_h1 = as.integer(unlist(column_order(h1)))


    probe_start = c(myProm$probeStart)
    enst_tss = unique(c(gene_data$transcript_start_site))
    enst_ids = character(length(enst_tss))
    enst_colors <- rainbow(length(enst_tss))

    probe_col <- rep("#FFFFFF", nrow(beta.data))
    enst_ids <- character(length(enst_tss))
    probe_pos <- numeric(length(enst_tss))
    closest_probe_ids <- character(length(unique(enst_tss)))

    # Main loop
    for (i in seq_along(enst_tss)) {
    tss_value <- enst_tss[i]
    closest_probe_index <- which.min(abs(probe_start - tss_value))
    closest_probe_start <- probe_start[closest_probe_index]
    probe_id <- myProm$probes[closest_probe_index]
    
    probe_pos[i] <- probe_id
    closest_probe_ids[i] <- probe_id
    
    enst_ids_pos <- which(gene_data$transcript_start_site == tss_value)[1]
    
    if (length(enst_ids_pos) > 1) {
        enst_ids[i] <- gene_data$transcript_id[enst_ids_pos[1]]
    } else {
        enst_ids[i] <- gene_data$transcript_id[enst_ids_pos]
    }
    if (is.na(probe_id)) {
        warning(paste("No matching probe_id found in beta data for", probe_id))
    } else{
        
        enst_index <- which(enst_ids == gene_data$transcript_id[enst_ids_pos])
    # cat("enst_pos: ", enst_index, "\n")
        probe_col[closest_probe_index] <- rainbow(length(enst_tss))[enst_index]
        names(probe_col)[closest_probe_index] <- enst_ids[i]
        }
    }

    # Assign non-empty names to probe_col
    probe_col_names <- names(probe_col)

    for (i in seq_along(probe_col)) {
    if (is.na(probe_col_names[i])) {
        probe_col_names[i] <- "NA"
    }
    }

    names(probe_col) <- probe_col_names

    probe_col_names <- ifelse(is.na(names(probe_col)), "NA", names(probe_col))

    # Convert IDs to characters
    enst <- as.character(sub(".*0{5}", "", probe_col_names))

    col_vector <- factor(enst, levels = unique(enst ))

    # Create and draw the heatmap
    h2 <- Heatmap(beta.data,
                col = my_palette, 
                name = "Methylation",
                cluster_rows = FALSE,
                column_order = col_order_h1,
                clustering_method_columns = "ward.D2",
                column_title = paste("Gene:", my.gene, "promoter: upS=",my.upS,"dnS=",my.dnS),
                show_heatmap_legend = TRUE,
                column_title_gp = gpar(fontsize = 10),
                column_names_gp = gpar(fontsize = 0.5),
                row_names_gp = gpar(fontsize = 5),
                width = 14,
                height = 10,
                left_annotation = HeatmapAnnotation(df = data.frame(enst), 
                                                    col = list(probe_col = col_vector), na_col = "white",
                                                    which = "row", show_annotation_name = FALSE, show_legend = TRUE))

    col_order_h2 = column_order(h2)


    # Check if column orders are identical 
    stopifnot(identical((col_order_h1), (col_order_h2)))

    # Arrange the heatmaps in a grid
    grob_h1 = ggplotify::as.grob(h1)
    grob_h2 = ggplotify::as.grob(h2)

    pdf(paste("D:/BINP37/DNA_repair_new/",my.gene,".pdf"))
    # Arrange the grobs in a grid
    grid.arrange(grob_h1, grob_h2, ncol = 1)
    dev.off()
}

