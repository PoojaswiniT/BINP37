###

# Description - the following code generates box plots for the FPKM values of FGFR4 gene
# and the methylation status of two key CpGs

library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)

# load the data
load("/Users/pooja/Documents/BINP37/Annotations_matched_EPIC_ERpHER2n.RData") # the annotation tracks provided
load("variance_data.Rdata") 
load("Documents/BINP37/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData")

#Load the function to split sample names
source("Documents/BINP37/sample_name_split.R")

#### Subselect FU=TRUE & ERpHER2n ##
ii=which( (pam50.frame$ClinGroup=="ERpHER2nLNn") | (pam50.frame$ClinGroup=="ERpHER2nLNp") )
pam50.frame=pam50.frame[ii,]

ii=which( pam50.frame$Follow.up.cohort  )
pam50.frame=pam50.frame[ii,]


### Load the BASE-processed EPIC data, edit names, match to WGS ###
load("Documents/BINP37/scanb_base.RData")
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

#Load the promoter annoataion object for EPIC data
load("Documents/BINP37/getPromoter.RData")

#Load the function to get find the probes
load("Documents/BINP37/get_cg.Rdata")

# Get the gene expression levels
gene_data = variance_data[variance_data$gene_symbol == "FGFR4", ]

# Remove the first few columns
gene_count_only = gene_data[,-c(1:11)]

gene_matrix = as.matrix(gene_count_only)

colnames(gene_matrix) = sapply(sample_name_split(colnames(gene_matrix)), `[`, 1)

rownames(gene_matrix) = c(gene_data$transcript_id)

# Summarise all the transcript expression to gene expression
gene_expression = colSums(gene_matrix, na.rm = TRUE)

# Create a data frame to store the gene expression data for each sample
gene_expression_df = data.frame(
  sample = names(gene_expression),
  gene_expression = gene_expression
)

#enst = gene_matrix["ENST00000502906",]

# Get the methylation data
my.upS = 10000
my.dnS = 2000
myProm = as.data.frame(get_cg(gene_data, my.upS = 10000, my.dnS = 2000))
myProm = myProm[!duplicated(myProm), ]

ii=match(myProm$probes,rownames(beta.matrix))
ii.i=which(!is.na(ii))
beta.data=beta.matrix[ii[ii.i],]

myProm$beta=beta.data

rownames(beta.data) = myProm$probes

beta.data = as.data.frame(t(beta.data))

beta.data = beta.data %>%
  select(cg00618323, cg20898288)

analysis_df <- data.frame(
  expression = gene_expression,
  cg00618323 = beta.data$cg00618323,
  cg20898288 = beta.data$cg20898288
)

# Convert to long format for plotting
analysis_long <- analysis_df %>%
  pivot_longer(
    cols = c(cg00618323, cg20898288),
    names_to = "probe",
    values_to = "methylation"
  )

analysis_long <- analysis_long %>%
  mutate(methylation_status = ifelse(methylation <= 0.4, "hypomethylated", "hypermethylated")) %>%
  mutate(methylation_status = factor(methylation_status, levels = c("hypomethylated", "hypermethylated")))

summary_data <- analysis_long %>%
  group_by(methylation_status) %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE),          # Mean expression
    se_expression = sd(expression, na.rm = TRUE) / sqrt(n())   # Standard error
  )

count_data <- analysis_long %>%
  group_by(methylation_status) %>%
  summarise(count = n())

# Plot the boxplot
ggplot(analysis_long, aes(x = methylation_status, y = expression)) +
  geom_boxplot(aes(fill = methylation_status)) +
  labs(x = "shore CpG clusters", y = "FPKM of FGFR4", title = "FGFR4 expression and methylation status") +
  theme_classic() +
  scale_fill_manual(values = c("hypomethylated" = "blue", "hypermethylated" = "yellow")) +
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        axis.ticks = element_line(color = "black"),  
        axis.line = element_line(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA)) + 
   geom_text(data = count_data, aes(x = methylation_status, y = max(analysis_long$expression) + 1, label = paste("n=",count)), 
            size = 5, color = "black", position = position_dodge(width = 0.75))
ggsave("Documents/BINP37/FGFR4_methylation.png")
