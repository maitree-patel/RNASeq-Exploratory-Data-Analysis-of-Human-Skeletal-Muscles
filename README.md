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
DESeq2 plots the log2FC for mean normalized counts through the plotMA() function, providing a distribution of the estimated log2FC across the genes giving a robust overview. Plot was visulaized using the more informative shrunken results.
![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/49c5f935-8819-4aca-94a1-7eba0de48930)
The plot shows the significant differentially expressed counts in blue.

### Counts Plots
Particular genes can be visualized by ploting counts for each comparison group, here, male and female. The gene with the highest log2FC (most differentially expressed) value was visualized.
![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/f3512c72-107c-4b7e-a480-fa94c0f4e53a)
The plot shows gene DEAD-box helicase 3 Y-linked (Ensembl: ENSG00000067048 and gene name: DDX3Y) which is most differentially expressed gene. This gene being Y-linked shows a significantly high expression level in males as would be expected. 


