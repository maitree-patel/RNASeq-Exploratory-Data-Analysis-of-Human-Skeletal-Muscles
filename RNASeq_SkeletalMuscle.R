#data: resting skeletal muscle biopies from untrained individuals (baseline transcriptome)
#question: to compare gene expression in skeletal muscles between male and females
##we know physiological difference exists

#loading the data 
#source: https://www.ebi.ac.uk/gxa/experiments/E-GEOD-58608/Downloads (EMBL expression atlas)
load("/Users/maitreepatel/Downloads/E-GEOD-58608-atlasExperimentSummary.Rdata")

#retrieving the count matrix
muscle_data <- experiment_summary@listData$rnaseq@assays$data@listData$counts
dim(muscle_data)

#retrieving the metadata
metadata.muscle <- experiment_summary@listData$rnaseq@colData
dim(metadata.muscle)

#checking to see if the number of samples are the same in the count matrix and the metadata
colnames(muscle_data) == rownames(metadata.muscle)

View(muscle_data)
View(as.data.frame(metadata.muscle))

library(DESeq2)

#creating a DESeq object
dds.obj <- DESeqDataSetFromMatrix(countData = muscle_data,
                                  colData = metadata.muscle,
                                  design = ~sex) #setting up our experimental design to compare between males and females

#performing DESeq analysis
dds.muscle <- DESeq(dds.obj)

#looking the dispersion estimate
plotDispEsts(dds.muscle)

#extracting result table from the DESeq2 analysis
res.muscle <- results(dds.muscle, alpha = 0.05)
summary(res.muscle)
View(as.data.frame(res.muscle))

#plotting log2 Fold Change versus normalized mean counts for each gene
plotMA(res.muscle)

#perfoming shrinking
resultsNames(dds.muscle)

lfc <- lfcShrink(dds = dds.muscle,
                 coef = 2, #has to be same as the design
                 type = "normal")

plotMA(lfc,
       ylim = c(-3,3))

#plotting the count of reads for the expression between male and females
plotCounts(dds.muscle, 
           gene=which.min(res.muscle$padj), 
           intgroup="sex")
#gene encodes for ribosomal protein S4 Y-linked 1 (RPS4Y1), a component of 40s ribosomal subunit. 
#This protein is encoded by this gene along with ribosomal protein S4, X-linked (RPS4X)
#it is characterised as cancer enhanced, linked to prostatde cancer. expected to see a higher expression level in males


RPS4X <- which(rownames(res.muscle) == "ENSG00000198034")
plotCounts(dds.muscle, 
           gene=RPS4X, 
           intgroup="sex")

#visualizing a cluster plot
sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)

library(RColorBrewer)

rownames(sampleDistMatrix) <- paste(colnames(vst), vst$sex, sep = "-")
colnames(sampleDistMatrix) <- paste(colnames(vst), vst$sex, sep = "-")
map_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

clustering_plot <- heatmap(sampleDistMatrix,
                           name = "Variance",
                           col = map_colors,
                           #row distance method
                           clustering_distance_rows = "pearson",
                           #column distance method
                           clustering_distance_columns = "pearson",
                           show_row_names = TRUE,
                           show_column_names = TRUE)

#heatmap of expression levels for the top 20 most expressed genes
#transforming data for better visualisation
vst <- varianceStabilizingTransformation(dds.muscle,
                                         blind = FALSE)

library(ComplexHeatmap)

genestokeep <- order(rowVars(assay(vst)), decreasing = T)[1:20]
heatmap_data <- as.data.frame(assay(vst)[genestokeep, ])

#getting the gene symbols
library(org.Hs.eg.db)

columns(org.Hs.eg.db)
heatmap_data$GENEID <- mapIds(org.Hs.eg.db,
                              keys = rownames(heatmap_data),
                              column = "GENENAME",
                              keytype = "ENSEMBL",
                              multiVals = "first")
heatmap_data$GENEID
#NAs found
heatmap_data$GENEID[12] <- "MTRNR2L8"
heatmap_data$GENEID[18] <- "TXLNGY"

#setting rownames to geneIDs
rownames(heatmap_data) <- heatmap_data$GENEID


#dropping the GENEID column before plotting
heatmap_data <- subset(heatmap_data,
                       select = -c(GENEID))    

#scaling to between -1 to 1 or centred around 0
heatmap_data2 = t(apply(heatmap_data, 1, function(x) {scale(x)}))
colnames(heatmap_data2) = colnames(heatmap_data)

# Create a HeatmapAnnotation object to specify the colours for conditions and replicates.
sex_colors = c("maroon", "dodgerblue")
names(sex_colors) = unique(metadata.muscle[,"sex"])

heatmap_anno = HeatmapAnnotation(Sex = as.matrix(colData(dds.muscle)[,c("sex")]), 
                                 col = list(Sex = sex_colors))

library(colorRamp2)
library(tidyverse)

expression_heatmap <- Heatmap(heatmap_data2,
                              name = "Expression",
                              top_annotation = heatmap_anno,
                              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                              cluster_rows = TRUE,
                              cluster_columns = FALSE,
                              #row distance method
                              clustering_distance_rows = "spearman",
                              #column distance method
                              clustering_distance_columns = "spearman",
                              #row cluster method
                              clustering_method_rows = "centroid",
                              #column cluster method
                              clustering_method_columns = "centroid",
                              row_names_side = "left",
                              show_row_names = TRUE,
                              show_column_names = TRUE,
                              column_names_side = "bottom",
                              #position of row dendrogam
                              row_dend_side = "left",
                              #position of column dendrogram
                              column_dend_side = "top",
                              #width of row dendrogram
                              row_dend_width = unit(1, "cm"),
                              #height of column dendrogram
                              column_dend_height = unit(1, "cm"))
expression_heatmap
#a distinct expression pattern between the genes that are y-linked s

#pca plot
pcaData <- plotPCA(vst,
                   intgroup = "sex",
                   returnData = TRUE)

percentVar <- round(100*attr(pcaData, "percentVar"))

pcaData %>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = sex)) +
  geom_point(size = 3) +
  labs(color = "sex") +
  xlab(paste0("PC1:", percentVar[1], "% variance")) +
  ylab(paste0("PC2:", percentVar[2], "% variance")) +
  theme_gray() +
  scale_color_manual(values = c("maroon", "dodgerblue"))

#volcano plot
res.df <- as.data.frame(res.muscle)
res.df <- mutate(res.df,
                 sig=ifelse(res.df$padj<0.1, "sig", "not_sig"))
View(res.df)  

res.df[which(abs(res.df$log2FoldChange)<1.0), "sig"] = "not_sig"

res.df %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_point(aes(col=sig)) +
  scale_colour_manual(values = c("red", "black"))

#volcano plot 2
res.df2 <- as.data.frame(res.muscle) 
res.df2$diffex <- "NO"
res.df2$diffex[res.df2$log2FoldChange>0.1 & res.df2$padj<0.05] <- "Upregulated"
res.df2$diffex[res.df2$log2FoldChange<0.1 & res.df2$padj<0.05] <- "Downregulated"

res.df2 %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_point(aes(col=diffex)) +
  scale_colour_manual(values = c("red", "black", "dodgerblue"))

#volcano plot 3
res.df3 <- as.data.frame(res.muscle)
res.df3$diffex <- "NO"
res.df3$diffex[res.df2$log2FoldChange>0.1] <- "Upregulated"
res.df3$diffex[res.df2$log2FoldChange<0.1] <- "Downregulated"
res.df3$sex <- "null"
res.df3[which(t(as.data.frame(metadata.muscle))[7,]=="male"), "sex"] <- "male"
res.df3[which(t(as.data.frame(metadata.muscle))[7,]=="female"), "sex"] <- "female"

res.df3 %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_point(aes(col=sex),
             alpha = 0.3) +
  scale_colour_manual(values = c("red", "black", "dodgerblue")) +
  ylim(0,25)

#uwsing enhancedcolcano package
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(lfc,
                lab = rownames(res.muscle),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,80),
                labSize = 4.0,
                FCcutoff = 0.5,
                col=c('black', 'green', 'blue', 'red3'),
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                              'p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0) 
#source: https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
