library(tidyverse)
library(ggpubr)
library(DESeq2)
library(ggplot2)
library(VennDiagram)
library(writexl)
library(pheatmap)

#R version 4.1.2 (2021-11-01)
#RStudio version 2021.09.0+351 
#tidyverse version 1.3.1 
#DESeq2 version 1.34.0
#BiocManager 1.30.16, R 4.1.1
#apeglm version 1.14.0
#ggplot2 version 3.3.5
#Writexl version 1.4.0
#pheatmap version 1.0.12
#ggpubr version 0.4.0

# Make a directory for plots: 'Plots'

#Produce PCA plot with DAL
#Import read counts
All <- read_tsv("./All.tsv",col_names = TRUE)

#Generate a matrix with the first row and first column as column/row names
All <- All %>% remove_rownames %>% column_to_rownames(var="gene")
cts <- as.matrix(All[ , ]) 

#Import sample details
coldata <- read_tsv("./Sample_annotation.tsv")
#Specify that the first column are the rownames
coldata <- data.frame(coldata, row.names = 1)

#Constructing a DESeq Dataset with TE1 as reference level
dds_TE1 <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata, 
  design= ~ Domain)
#Set TE1 as reference level
dds_TE1$Domain <- relevel(dds_TE1$Domain, ref="TE1")
#Differential expression analysis with TE1 as reference
dds_TE1 <- DESeq(dds_TE1)
resultsNames(dds_TE1)

#Plot how related the samples are to eachother
vstcounts <- vst(dds_TE1, blind=TRUE)
PCAplot <- plotPCA(vstcounts, intgroup=c("Domain"), returnData=TRUE)
percentVar <- round(100 * attr(PCAplot, "percentVar"))
PCAplot2 <- ggplot(PCAplot, aes(PC1, PC2, shape=Domain, color=Domain)) +
  geom_point(size=2.5, alpha = 1) +
  scale_shape_manual(values=c(6, 4, 2, 5)) +  
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +  
  theme(legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(legend.position = c(0.25, 0.4)) + 
  theme(legend.title = element_text(size=8))+ 
  theme(legend.text= element_text(size=8)) 

PCAplot2
ggsave('Figure 2A Endosperm domains PCA plot.pdf', width=7.78, height=2.67, plot=PCAplot2, device = "pdf", path="./Plots/", dpi=600)

#Produce MAplots without DAL
#Import read counts
All2 <- read_tsv("./All.tsv",col_names = TRUE)
All2 <- subset(All2, select=-c(DAL1_5_2pos, DAL1_5_4pos))

#Generate a matrix with the first row and first column as column/row names
All2 <- All2 %>% remove_rownames %>% column_to_rownames(var="gene")
cts <- as.matrix(All2[ , ]) 

#Import sample details
coldata2 <- read_tsv("./Sample_annotation_2.tsv")
#Specify that the first column are the rownames
coldata2 <- data.frame(coldata2, row.names = 1)

#Constructing a DESeq Dataset with TE1 as reference level
dds_TE1 <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata2, 
  design= ~ Domain)
#Set TE1 as reference level
dds_TE1$Domain <- relevel(dds_TE1$Domain, ref="TE1")
#Differential expression analysis with TE1 as reference
dds_TE1 <- DESeq(dds_TE1)
resultsNames(dds_TE1)

#Differential gene expression
EE_vs_TE1 <- lfcShrink(dds_TE1, coef=2, type="apeglm")
ESR_vs_TE1 <- lfcShrink(dds_TE1, coef=3, type="apeglm")

summary(EE_vs_TE1)
#MAplots
#Export MAplot with only a line at log2FC = 1
EE_vs_TE1_MAplot <- ggmaplot(
  EE_vs_TE1, 
  top=0,
  fdr = 0.05, fc = 1, size = 0.05, 
  palette = c("#B31B21", "#1465AC", "darkgray"), 
  xlab = "Log2 Mean Expression",  ylab="Informative read Log2 FC EE vs TE1") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(legend.position="top") +
  geom_hline(yintercept=c(-1,1), linetype="dotted") +
  annotate("text", x=0.2, y=16, label="ns: 8774", size=4) +
  scale_x_continuous(breaks = seq(0, 18, by= 2), limits=c(0, 18)) +
  scale_y_continuous(breaks = seq(-20, 20, by = 4), limits=c(-22, 22)) 
  
EE_vs_TE1_MAplot
ggsave('SFigure 2B EE_vs_TE1.pdf', width=9, height=9, plot=EE_vs_TE1_MAplot, device = "pdf", path="./Plots/MAplots/", dpi=1200)

ESR_vs_TE1_MAplot <- ggmaplot(
  ESR_vs_TE1, 
  top=0,
  fdr = 0.05, fc = 1, size = 0.05, 
  palette = c("#B31B21", "#1465AC", "darkgray"), 
  xlab = "Log2 Mean Expression",  ylab="Informative read Log2 FC ESR vs TE1") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(legend.position="top") +
  geom_hline(yintercept=c(-1,1), linetype="dotted") +
  annotate("text", x=0.2, y=11, label="ns: 10427", size=4) +
  scale_x_continuous(breaks = seq(0, 18, by= 2), limits=c(0, 18)) +
  scale_y_continuous(breaks = seq(-12, 12, by = 4), limits=c(-12, 12)) 

ESR_vs_TE1_MAplot
ggsave('Figure 2B ESR_vs_TE1.pdf', width=9, height=6, plot=ESR_vs_TE1_MAplot, device = "pdf", path="./Plots/MAplots/", dpi=1200)

#EE_vs_ESR
#Constructing a DESeq Dataset with ESR as reference level
dds_ESR <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata2, 
  design= ~ Domain)
#Set TE1 as reference level
dds_ESR$Domain <- relevel(dds_ESR$Domain, ref="ESR")
#Differential expression analysis with ESR as reference
dds_ESR <- DESeq(dds_ESR)
resultsNames(dds_ESR)

#Differential gene expression
EE_vs_ESR <- lfcShrink(dds_ESR, coef=2, type="apeglm")

EE_vs_ESR_MAplot <- ggmaplot(
  EE_vs_ESR, 
  top=0,
  fdr = 0.05, fc = 1, size = 0.05, 
  palette = c("#B31B21", "#1465AC", "darkgray"), 
  xlab = "Log2 Mean Expression",  ylab="Informative read Log2 FC EE vs ESR") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(legend.position="top") +
  geom_hline(yintercept=c(-1,1), linetype="dotted") +
  annotate("text", x=0.2, y=16, label="ns: 7035", size=4) +
  scale_x_continuous(breaks = seq(0, 18, by= 2), limits=c(0, 18)) +
  scale_y_continuous(breaks = seq(-20, 20, by = 4), limits=c(-22, 22)) 

EE_vs_ESR_MAplot
ggsave('SFigure 2A EE_vs_ESR.pdf', width=9, height=9, plot=EE_vs_ESR_MAplot, device = "pdf", path="./Plots/MAplots/", dpi=1200)

ESR_vs_TE1_reads <- All2 %>% select(ESR2pos, ESR3pos, ESR4pos, TE1_2pos, TE1_3pos, TE1_4pos)
ESR_vs_TE1_reads[is.na(ESR_vs_TE1_reads)] <- 0
ESR_vs_TE1_reads <- tibble::rownames_to_column(ESR_vs_TE1_reads, "Gene")
#Import the full DEseq per comparison
ESR_vs_TE1 <- as.data.frame(ESR_vs_TE1)
ESR_vs_TE1 <- tibble::rownames_to_column(ESR_vs_TE1, "Gene")
#Pull out log2FC and padj columns
ESR_vs_TE1 <- ESR_vs_TE1 %>% select("Gene", log2FoldChange, padj)
#Combine reads with log2FC and padj columns
ESR_vs_TE1_reads <- full_join(ESR_vs_TE1_reads, ESR_vs_TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 
ESR_vs_TE1_reads <- ESR_vs_TE1_reads[order(ESR_vs_TE1_reads$log2FoldChange),]

EE_vs_TE1_reads <- All2 %>% select(EE1pos, EE3pos, EE4pos, TE1_2pos, TE1_3pos, TE1_4pos)
EE_vs_TE1_reads[is.na(EE_vs_TE1_reads)] <- 0
EE_vs_TE1_reads <- tibble::rownames_to_column(EE_vs_TE1_reads, "Gene")
#Import the full DEseq per comparison
EE_vs_TE1 <- as.data.frame(EE_vs_TE1)
EE_vs_TE1 <- tibble::rownames_to_column(EE_vs_TE1, "Gene")
#Pull out log2FC and padj columns
EE_vs_TE1 <- EE_vs_TE1 %>% select("Gene", log2FoldChange, padj)
#Combine reads with log2FC and padj columns
EE_vs_TE1_reads <- full_join(EE_vs_TE1_reads, EE_vs_TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 
EE_vs_TE1_reads <- EE_vs_TE1_reads[order(EE_vs_TE1_reads$log2FoldChange),]

EE_vs_ESR_reads <- All2 %>% select(EE1pos, EE3pos, EE4pos, ESR2pos, ESR3pos, ESR4pos)
EE_vs_ESR_reads[is.na(EE_vs_ESR_reads)] <- 0
EE_vs_ESR_reads <- tibble::rownames_to_column(EE_vs_ESR_reads, "Gene")
#Import the full DEseq per comparison
EE_vs_ESR <- as.data.frame(EE_vs_ESR)
EE_vs_ESR <- tibble::rownames_to_column(EE_vs_ESR, "Gene")
#Pull out log2FC and padj columns
EE_vs_ESR <- EE_vs_ESR %>% select("Gene", log2FoldChange, padj)
#Combine reads with log2FC and padj columns
EE_vs_ESR_reads <- full_join(EE_vs_ESR_reads, EE_vs_ESR, by="Gene", copy=FALSE,suffix=c("1","2")) 
EE_vs_ESR_reads <- EE_vs_ESR_reads[order(EE_vs_ESR_reads$log2FoldChange),]

# Make a directory for data output: 'Output'
write_xlsx(list("PCA plot input"=All, "MAplot EE vs TE1"=EE_vs_TE1_reads, "MAplot ESR vs TE1"=ESR_vs_TE1_reads, "MAplot EE vs ESR"=EE_vs_ESR_reads), "./Output/SData2 Domains_Differential_Gene_Expression.xlsx")
