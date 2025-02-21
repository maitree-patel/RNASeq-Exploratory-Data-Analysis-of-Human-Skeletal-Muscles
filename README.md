# RNASeq Exploratory Data Analysis of Human Skeletal Muscles
RNASeq data of human skeletal muscles was to analyzed to asses the differences between males and females gene expression

## Objective
A Publication titled "The human skeletal muscle transcriptome: sex differences, alternative splicing, and tissue homogeneity assessed with RNA sequencing" extensively explored the baseline transcriptome of the human skeletal muscles, comapring differences in males in females. The main objective of this project was to load, tranform, visualize and interpret RNASeq data made available through this study. 

## Data and Methodology
RNASeq data was accessed through the EMBL expression atlas. The data is also available at the Gene Expression Omnibus (GEO) under the accession number GSE58387 and GSE58608. The data consisted of 26 total male and female samples taken from untrained skeletal muscles from the leg and a total of 58735 gene were assesed.

DESeq2 package from R was used to perform differential gene expression analysis followed by downstream manipulation and visualizations to exploratorily analyze data, particularly for visualization. DeSeq2 uses log2FoldChange to estimate differential expression. tidyverse, ggplot2, CompleHeatmap and EnhancedVolcano packages were additionally used to manipulate and visualize data. The code has been made available as a part of this repository.

## Exploring Results
Shrinkage of effect size (LFC estimate) was performed on the results from the DESeq analysis to account for the high variation in expression levels of lowly expressed gene. With a default p-value of 0.01, a total of 636 out of 32588 (2%) genes were upregulated and 472 (1.4%) were downregulated. 

### Dispersion Estimates Plot
For the dispersion estimates for each gene as measured by DESeq, a curve is fit through it to estimate dispersion levels based on gene expression. This is plotted using the plotDispEsts() function.

![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/52050aa3-a4d9-4404-9a9c-23cafa254000)

DESeq2 assumes that genes with similar expression levels have similar dispersion, and this used to predict the variation of genes around the mean or dispersion. The red line is the curve that estimates the dispersion for each gene based on the fitted data. The plot follows the expected curve, i.e. the dispersion decreses as the gene expression increases (mean normalized count).

### MA-Plot
DESeq2 plots the log2FC for mean normalized counts through the plotMA() function, providing a distribution of the estimated log2FC across the genes giving a robust overview. Plot was visulaized using the more informative shrunken results.

![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/49c5f935-8819-4aca-94a1-7eba0de48930)

The plot shows the significant differentially expressed counts in blue.

### Counts Plots
Particular genes can be visualized by ploting counts for each comparison group, here, male and female. The gene with the highest log2FC (most differentially expressed) value was visualized.

![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/f3512c72-107c-4b7e-a480-fa94c0f4e53a)

The plot shows gene DEAD-box helicase 3 Y-linked (Ensembl: ENSG00000067048 and gene name: DDX3Y) which is most differentially expressed gene. This gene being Y-linked shows a significantly high expression level in males as would be expected. 

The gene most significantly expressed (lowest p adjusted value was also visualized.

![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/f1e07cbb-fa36-4452-865f-55e30d49b6d2)

The plot shows another Y-linked gene called ribosomal protein S4 Y-linked 1, which encodes protein component of 40s ribosomal subunit. It is characterised as cancer enhanced gene, linked to prostatde cancer and therefore also expected to see a higher expression level in males as evident through the visualization as well.

Since females have two X-chromosomes (one silenced for dosage compensation) and males have a single X-chromosome along with theY chromosomel, both would have the same number of X-linked genes. This in evident through the below count plot that visualzes one of the X-linked genes, ribosomal protein S4 X-linked, which along with the ribosomal protein S4 Y-linked 1 is responsible for encoding the 40s ribosomal subunit.

![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/c8c7e73b-01c0-4f1c-ab68-140ded54bb92)

### Clustering Map
Similarity between the samples were assessed by measuring distances between them and visualizing a cluster heatmap.

![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/24a039cb-ca10-4c65-b140-6a33c7b0c080)

We can see that majorly the male samples cluster closer together with other male samples, same being true for the female samples.

### PCA Plot
The PCA below shows similar results as the clustering, i.e. distance between the samples, where we see a clear clustering of male samples together and away from the clustering of the female samples.

![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/b85f9352-391e-4e1f-b2f5-54dadba2674b)

### Count Matrix Heatmap
Heatmap showing expression levels of the top 20 genes were visualized using the Heatmap() function in the ComplexHeatmap package.

 ![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/cd7ac4a9-2e34-4610-b71c-4265a3a1d9c9)

- We can see a large cluster of genes starting from UTY until RPS4Y1 which are Y-linked genes seen to be highly upregulated in male samples with distinct differentiation from the female samples that have them downregulated.
- The MAOA gene shows upregulation in 5 out of the total 13 male samples and 3 out of the total 13 female sample. MAOA is the Monoamine oxidase A gene, also known as the "warrior gene" which has been linked to behavioural agression.
- We see the upregulation of the XIST, i.e. X inactive specific transcript gene in all the female samples and a downregulation in the male sample. This gene is involved in an early developmental process of the silencing of one of X chromosome in a pair in females. The gene is linked to the skewed X-linked inactivation disease that favours the inactivation of one chromosome over the other and has been shown to interact with BRCA1, i.e. Breast Cancer 1 Suseptibility protein gene.

### Volcano Plot
The volcano plot are essential visualizations that provide an overview of the expression of genes in the samples. It plots the adjusted p-values, a significance measure, on the Y-axis and the log2FC on the X-axis. The Volcano plot below shows the significantly down- and upregulated genes with the top 20 genes labeled.

<img width="926" alt="Screenshot 2023-12-23 at 9 25 21 PM" src="https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/0155f333-0696-4c6b-bd91-8ad17db05bb7">

The EnhancedVolcano package was used to exploratorily analyze visualizing the same volcano plot color coded with different conditions.

![image](https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles/assets/134908239/ac6f2795-6b61-429f-ba47-b8b34340dde4)

## References
1. https://faseb.onlinelibrary.wiley.com/doi/epdf/10.1096/fj.14-255000
2. https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts
3. https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
4. https://www.ebi.ac.uk/gxa/experiments/E-GEOD-58608/Experiment%20Design
5. https://www.pnas.org/doi/10.1073/pnas.0808376106
6. https://uclouvain-cbio.github.io/WSBIM2122/sec-rnaseq.html











