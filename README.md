# RNASeq Exploratory Data Analysis of Human Skeletal Muscles
RNASeq data of human skeletal muscles was to analyzed to asses the differences between males and females gene expression

## Objective
A Publication titled "The human skeletal muscle transcriptome: sex differences, alternative splicing, and tissue homogeneity assessed with RNA sequencing" extensively explored the baseline transcriptome of the human skeletal muscles, comapring differences in males in females. The main objective of this project was to load, tranform, visualize and interpret RNASeq data made available through this study. 

## Data and Methodology
RNASeq data was accessed through the EMBL expression atlas. The data is also available at the Gene Expression Omnibus (GEO) under the accession number GSE58387 and GSE58608. The data consisted of 26 total male and female samples taken from untrained skeletal muscles from the leg and a total of 58735 gene were assesed.

DESeq2 package from R was used to perform differential gene expression analysis followed by downstream manipulation and visualizations to exploratorily analyze data, particularly for visualization. DeSeq2 uses log2FoldChange to estimate differential expression. tidyverse, ggplot2, CompleHeatmap and EnhancedVolcano packages were additionally used to manipulate and visualize data. The code has been made available as a part of this repository.

## Exploring Results
Shrinkage of effect size (LFC estimate) was performed on the results from the DESeq analysis to account for the high variation in expression levels of lowly expressed gene. With a default p-value of 0.01, a total of 636 out of 32588 (2%) genes were upregulated and 472 (1.4%) were downregulated. 

### MA-Plot
![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/49c5f935-8819-4aca-94a1-7eba0de48930)
