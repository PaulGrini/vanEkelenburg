#Differential gene expression between Col-0 and Tsu-1 at 4 DAP (Early)
#Collect output data in '/Output'
#Collect output plots in '/Plots'

library(tidyverse)
library(ggpubr)
library("DESeq2")
library("ggplot2")
library(VennDiagram)
#R version 4.1.1 (2021-08-10)
#RStudio version 1.4.1106 
#tidyverse version 1.3.1 
#DESeq2 version 1.32.0
#BiocManager 1.30.16, R 4.1.0
#apeglm version 1.14.0
#ggplot2 version 3.3.4
#mvtnorm version 1.1-2

#Import reads
Early_All <- read_tsv("./Merged_read_counts/Early_all_reads.tsv",col_names = TRUE)

#Generate a matrix with the first row and first column as column/row names
Early_All <- Early_All %>% remove_rownames %>% column_to_rownames(var="gene")
cts_Early_All <- as.matrix(Early_All[ , ]) 

#Import sample details
coldata_col_vs_Tsu <- read_tsv("./Sample_annotation/Sample_annotation_ColTsu_early.tsv")
#Specify that the first column are the rownames
coldata_col_vs_Tsu <- data.frame(coldata_col_vs_Tsu, row.names = 1)

#Boxplot of read count distribution per sample
boxplot(log2(cts_Early_All + 1))

#Constructing a generic DESeq Dataset for Early
dds_EarlyColTsu <- DESeqDataSetFromMatrix(
  countData = cts_Early_All,
  colData = coldata_col_vs_Tsu, 
  design= ~ Ecotype)
#Differential expression analysis with Col as reference
dds_EarlyColTsu$Ecotype <- relevel(dds_EarlyColTsu$Ecotype, ref="Tsu")
dds_EarlyColTsu <- DESeq(dds_EarlyColTsu)
resultsNames(dds_EarlyColTsu)

EarlyColTsu_DESeq <- lfcShrink(dds_EarlyColTsu, coef=2, type="apeglm")

#Export DESeq2 data
Early_Col_vs_Tsu <- EarlyColTsu_DESeq[order(EarlyColTsu_DESeq$log2FoldChange),] 
write.csv(as.data.frame(Early_Col_vs_Tsu), file="./Output/Early_Col_vs_Tsu.csv")

summary(EarlyColTsu_DESeq, alpha = 0.05)

#Plot how related the samples are to eachother (PCAplot)
vstcounts <- vst(dds_EarlyColTsu, blind=TRUE)
pdf('./Plots/PCAplot_Early_Col_vs_Tsu.pdf')
plotPCA(vstcounts, intgroup=c("Ecotype"))
dev.off()

#Export MAplot
#Export MAplot with only a line at log2FC = 1
maplot_lfc1 <- ggmaplot(
  EarlyColTsu_DESeq, 
  top=0,
  fdr = 0.05, fc = 1, size = 0.05, 
  palette = c("#B31B21", "#1465AC", "darkgray"), 
  xlab = "Log2 Mean Expression",  ylab="Log2 FC") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(legend.position="top") +
  geom_hline(yintercept=c(-1,1), linetype="dotted") +
  annotate("text", x=0.2, y=1.5, label="log2FC = 1", size=2) +
  annotate("text", x=0.2, y=-0.5, label="log2FC = -1", size=2) + 
  scale_x_continuous(breaks = seq(0, 21, by= 2), limits=c(0, 21)) +
  scale_y_continuous(breaks = seq(-22, 22, by = 2), limits=c(-22, 22)) 

maplot_lfc1
ggsave('SFigure 4A Early_Col_vs_Tsu_MAplot.pdf', width=10.6414, height=9, plot=maplot_lfc1, device = "pdf", path="./Plots/", dpi=1200)

