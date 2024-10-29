# Visualising the correlation of Isoform expression with DNA methylation in breast cancer cells. 

This repository contains a set of R scripts designed to analyze and visualize the correlation between isoform expression and DNA methylation in breast cancer cells. The provided scripts allow for the filtering, mapping, and visualization of data, ultimately producing heatmaps that illustrate these correlations.

## Table of Contents

- [Overview](#overview)
- [Installation and Requirements](#installation)
- [Content](#content)
  - [1. Mapping](#1-mapping)
  - [2. Filtering](#2-filtering)
  - [3. Heatmap Generation](#3-heatmap-generation)
  - [4. Gene Set Enrichment Analysis](#4-gene-set-enrichment-analysis)
- [Additional Scripts](#additional-scripts)
- [Usage](#usage)
- [Acknowledgements](#acknowledgements)

## Overview

This project focuses on understanding the relationship between isoform expression and DNA methylation in breast cancer cells. By using the provided transcript count matrix, SCAN-B dataset and methylation data, the scripts filter relevant data, map it to specific genomic annotations, and visualize the correlations through heatmaps. The analysis includes filtering transcripts based on transcript variance, FPKM values, isoform expression levels and DNA methylation levels, which are key in understanding breast cancer subtypes and progression.

## Installation and requirements

The code was executed in R Studio (v4.3.3)

The following R packages are required to run the scripts:

- **ggplot2**
- **pheatmap**
- **GSEA**
- **dplyr**
- **tidyr**
- **biomaRt**
- **GenomicRanges**
- **gplots**
- **RColorBrewer**
- **grid**
- **gridExtra**
- **ComplexHeatmap**

You can install these packages using the following command in R:

```r
install.packages(c("ggplot2", "pheatmap", "GSEA", "dplyr", "tidyr", "biomaRt", "GenomicRanges", "gplots", "RColorBrewer", "RColorBrewer","grid", "gridExtra", "ComplexHeatmap"))
```

## Content

### 1. Mapping

The `1_mapping.R` script filters the final sample list based on a provided annotation file. It narrows down the large set of samples to 515 specific samples and maps them to their respective ENsG, ENST IDs, and other attributes.

### 2. Filtering

The `2_filtering.R` script performs arbitrary filtering steps to refine the dataset. Genes with less than two transcripts, transcript variance of less than 10, FPKM value less than 1, and percentage of samples with less than 10% variance are filtered out.

### 3. Heatmap Generation

The `3_heatmap.R` script generates heatmaps by gathering the beta matrix from the SCAN-B database and a provided short script. It plots DNA methylation levels and isoform expression levels, including the transcription start site (TSS).

### 4. Gene Set Enrichment Analysis

The `GSEA.R` script performs Gene Set Enrichment Analysis, allowing for a deeper understanding of the functional significance of the identified genes and transcripts.


## Additional Scripts

- **CpG Mapping:** Additional scripts map CpG IDs to the beta data.
- **Box Plots:** Generate box plots for individual genes correlating with breast cancer subtypes.
- **Functions:** Various functions are included to facilitate running the main scripts.


## Usage 
Each script can be run in a sequence, but be sure to save the R data objects whenever prompted, as this will make it easier to run the scripts and continue the analysis smoothly. Also, all the dependencies must be installed and the input data must be formatted correctly, additional comments are provided in the scripts.

## Acknowledgments

I would like to thank Johan Staaf (Division of Translational Cancer Research) and team for their guidance and support throughout this project.