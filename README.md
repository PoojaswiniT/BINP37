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
