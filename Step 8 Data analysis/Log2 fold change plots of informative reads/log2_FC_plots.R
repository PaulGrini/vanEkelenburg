library(tidyverse)
library(ggplot2)
library(VennDiagram)
library(writexl)
#R version 4.1.1 (2021-08-10)
#RStudio version 2021.09.0+351 
#tidyverse version 1.3.1
#ggplot2 version 3.3.5
#VennDiagram version 1.6.20
#Writexl version 1.4.0

#Set working directory
setwd("/Users/yurisv/Desktop/Manuscripts/Domain paper/Heterozygous reads/IRP")

#Read all tables
EE <- read.csv("../Total inf read and ecotype specific gene filter/EE_inf_reads_homozygous_reads_and_gene_description.csv", header=TRUE)
names(EE)[names(EE) == "Col1"] <- "Col1EE"
names(EE)[names(EE) == "Col2"] <- "Col2EE"
names(EE)[names(EE) == "Col3"] <- "Col3EE"
names(EE)[names(EE) == "Tsu1"] <- "Tsu1EE"
names(EE)[names(EE) == "Tsu2"] <- "Tsu2EE"
names(EE)[names(EE) == "Tsu3"] <- "Tsu3EE"
ESR <- read.csv("../Total inf read and ecotype specific gene filter/ESR_inf_reads_homozygous_reads_and_gene_description.csv", header=TRUE)
names(ESR)[names(ESR) == "Col1"] <- "Col1ESR"
names(ESR)[names(ESR) == "Col2"] <- "Col2ESR"
names(ESR)[names(ESR) == "Col3"] <- "Col3ESR"
names(ESR)[names(ESR) == "Tsu1"] <- "Tsu1ESR"
names(ESR)[names(ESR) == "Tsu2"] <- "Tsu2ESR"
names(ESR)[names(ESR) == "Tsu3"] <- "Tsu3ESR"
TE1 <- read.csv("../Total inf read and ecotype specific gene filter/TE1_inf_reads_homozygous_reads_and_gene_description.csv", header=TRUE)
names(TE1)[names(TE1) == "Col1"] <- "Col1TE1"
names(TE1)[names(TE1) == "Col2"] <- "Col2TE1"
names(TE1)[names(TE1) == "Col3"] <- "Col3TE1"
names(TE1)[names(TE1) == "Tsu1"] <- "Tsu1TE1"
names(TE1)[names(TE1) == "Tsu2"] <- "Tsu2TE1"
names(TE1)[names(TE1) == "Tsu3"] <- "Tsu3TE1"

#Adj.Pvalue < 0.05
EE_MEG <- EE[EE$Informative_read_log2FC >= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
ESR_MEG <- ESR[ESR$Informative_read_log2FC >= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
TE1_MEG <- TE1[TE1$Informative_read_log2FC >= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),]
#PEGs
#Adj.Pvalue < 0.05 
EE_PEG <- EE[EE$Informative_read_log2FC <= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
ESR_PEG <- ESR[ESR$Informative_read_log2FC <= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
TE1_PEG <- TE1[TE1$Informative_read_log2FC <= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),]

#read counts of TE1 specific MEGs/PEGs in ESR and EE 
TE1_MEGs_in_ESR <- subset(ESR, (Gene %in% TE1_MEG$Gene))
TE1_PEGs_in_ESR <- subset(ESR, (Gene %in% TE1_PEG$Gene))
TE1_MEGs_in_EE <- subset(EE, (Gene %in% TE1_MEG$Gene))
TE1_PEGs_in_EE <- subset(EE, (Gene %in% TE1_PEG$Gene))

#TE1 specific MEGs/PEGs that are not present in EE or ESR dataset
TE1_MEGs_not_in_ESR_dataset <- subset(TE1_MEG, !(Gene %in% ESR$Gene))
TE1_PEGs_not_in_ESR_dataset <- subset(TE1_PEG, !(Gene %in% ESR$Gene))
TE1_MEGs_not_in_EE_dataset <- subset(TE1_MEG, !(Gene %in% EE$Gene))
TE1_PEGs_not_in_EE_dataset <- subset(TE1_PEG, !(Gene %in% EE$Gene))

#read counts of EE specific MEGs/PEGs in ESR and TE1
EE_MEGs_in_TE1 <- subset(TE1, (Gene %in% EE_MEG$Gene))
EE_PEGs_in_TE1 <- subset(TE1, (Gene %in% EE_PEG$Gene))

#EE specific MEGs/PEGs that are not present in TE1 or ESR dataset
EE_MEGs_not_in_TE1_dataset <- subset(EE_MEG, !(Gene %in% TE1$Gene))
EE_PEGs_not_in_TE1_dataset <- subset(EE_PEG, !(Gene %in% TE1$Gene))

#read counts of ESR specific MEGs/PEGs in TE1 and EE 
ESR_MEGs_in_TE1 <- subset(TE1, (Gene %in% ESR_MEG$Gene))
ESR_PEGs_in_TE1 <- subset(TE1, (Gene %in% ESR_PEG$Gene))

#ESR specific MEGs/PEGs that are not present in EE or TE1 dataset
ESR_MEGs_not_in_TE1_dataset <- subset(ESR_MEG, !(Gene %in% TE1$Gene))
ESR_PEGs_not_in_TE1_dataset <- subset(ESR_PEG, !(Gene %in% TE1$Gene))

#List of MEGs and PEGs for each domain which shows foldchanges in the other domains

#All ESR MEGs with log2foldchanges in TE1
ESR_MEGs_log2FC <- ESR_MEG[,c(1:8,10)]
ESR_MEGs_log2FC <- ESR_MEGs_log2FC[order(-ESR_MEGs_log2FC$Informative_read_log2FC),]
ESR_MEGs_TE1_log2FC <- ESR_MEGs_in_TE1[,c(1:8,10)]
names(ESR_MEGs_log2FC)[names(ESR_MEGs_log2FC) == "Informative_read_log2FC"] <- "log2FC_ESR"
names(ESR_MEGs_TE1_log2FC)[names(ESR_MEGs_TE1_log2FC) == "Informative_read_log2FC"] <- "log2FC_TE1"
names(ESR_MEGs_log2FC)[names(ESR_MEGs_log2FC) == "Adj.Pvalue"] <- "Adj.P_ESR"
names(ESR_MEGs_TE1_log2FC)[names(ESR_MEGs_TE1_log2FC) == "Adj.Pvalue"] <- "Adj.P_TE1"
ESR_MEGs_in_TE1 <- full_join(ESR_MEGs_log2FC, ESR_MEGs_TE1_log2FC, by="Gene", copy=FALSE,suffix=c("1","2")) 
ESR_MEGs_in_TE1$log2FC_TE1 <- as.character(ESR_MEGs_in_TE1$log2FC_TE1)
ESR_MEGs_in_TE1$log2FC_TE1[is.na(ESR_MEGs_in_TE1$log2FC_TE1)] <- "Not present"

#All ESR PEGs with log2foldchanges in TE1
ESR_PEGs_log2FC <- ESR_PEG[,c(1:8,10)]
ESR_PEGs_log2FC <- ESR_PEGs_log2FC[order(-ESR_PEGs_log2FC$Informative_read_log2FC),]
ESR_PEGs_TE1_log2FC <- ESR_PEGs_in_TE1[,c(1:8,10)]
names(ESR_PEGs_log2FC)[names(ESR_PEGs_log2FC) == "Informative_read_log2FC"] <- "log2FC_ESR"
names(ESR_PEGs_TE1_log2FC)[names(ESR_PEGs_TE1_log2FC) == "Informative_read_log2FC"] <- "log2FC_TE1"
names(ESR_PEGs_log2FC)[names(ESR_PEGs_log2FC) == "Adj.Pvalue"] <- "Adj.P_ESR"
names(ESR_PEGs_TE1_log2FC)[names(ESR_PEGs_TE1_log2FC) == "Adj.Pvalue"] <- "Adj.P_TE1"
ESR_PEGs_in_TE1 <- full_join(ESR_PEGs_log2FC, ESR_PEGs_TE1_log2FC, by="Gene", copy=FALSE,suffix=c("1","2"))
ESR_PEGs_in_TE1$log2FC_TE1 <- as.character(ESR_PEGs_in_TE1$log2FC_TE1)
ESR_PEGs_in_TE1$log2FC_TE1[is.na(ESR_PEGs_in_TE1$log2FC_TE1)] <- "Not present"

#All ESR MEGs and PEGs with the log2foldchanges for each gene in the other domains
#If gene is not present, this is  highlighted
ESR_all <- rbind(ESR_MEGs_in_TE1, ESR_PEGs_in_TE1)

#EE specific MEGs and PEGs
EE_MEGs_log2FC <- EE_MEG[,c(1:8,10)]
EE_MEGs_log2FC <- EE_MEGs_log2FC[order(-EE_MEGs_log2FC$Informative_read_log2FC),]
EE_MEGs_TE1_log2FC <- EE_MEGs_in_TE1[,c(1:8,10)]
names(EE_MEGs_log2FC)[names(EE_MEGs_log2FC) == "Informative_read_log2FC"] <- "log2FC_EE"
names(EE_MEGs_TE1_log2FC)[names(EE_MEGs_TE1_log2FC) == "Informative_read_log2FC"] <- "log2FC_TE1"
names(EE_MEGs_log2FC)[names(EE_MEGs_log2FC) == "Adj.Pvalue"] <- "Adj.P_EE"
names(EE_MEGs_TE1_log2FC)[names(EE_MEGs_TE1_log2FC) == "Adj.Pvalue"] <- "Adj.P_TE1"
EE_MEGs_in_TE1 <- full_join(EE_MEGs_log2FC, EE_MEGs_TE1_log2FC, by="Gene", copy=FALSE,suffix=c("1","2"))
EE_MEGs_in_TE1$log2FC_TE1 <- as.character(EE_MEGs_in_TE1$log2FC_TE1)
EE_MEGs_in_TE1$log2FC_TE1[is.na(EE_MEGs_in_TE1$log2FC_TE1)] <- "Not present"

EE_PEGs_log2FC <- EE_PEG[,c(1:8,10)]
EE_PEGs_log2FC <- EE_PEGs_log2FC[order(-EE_PEGs_log2FC$Informative_read_log2FC),]
EE_PEGs_TE1_log2FC <- EE_PEGs_in_TE1[,c(1:8,10)]
names(EE_PEGs_log2FC)[names(EE_PEGs_log2FC) == "Informative_read_log2FC"] <- "log2FC_EE"
names(EE_PEGs_TE1_log2FC)[names(EE_PEGs_TE1_log2FC) == "Informative_read_log2FC"] <- "log2FC_TE1"
names(EE_PEGs_log2FC)[names(EE_PEGs_log2FC) == "Adj.Pvalue"] <- "Adj.P_EE"
names(EE_PEGs_TE1_log2FC)[names(EE_PEGs_TE1_log2FC) == "Adj.Pvalue"] <- "Adj.P_TE1"
EE_PEGs_in_TE1 <- full_join(EE_PEGs_log2FC, EE_PEGs_TE1_log2FC, by="Gene", copy=FALSE,suffix=c("1","2")) 
EE_PEGs_in_TE1$log2FC_TE1 <- as.character(EE_PEGs_in_TE1$log2FC_TE1)
EE_PEGs_in_TE1$log2FC_TE1[is.na(EE_PEGs_in_TE1$log2FC_TE1)] <- "Not present"

EE_all <- rbind(EE_MEGs_in_TE1, EE_PEGs_in_TE1)

#TE1 specific MEGs and PEGs in EE
TE1_MEGs_log2FC <- TE1_MEG[,c(1:8,10)]
TE1_MEGs_log2FC <- TE1_MEGs_log2FC[order(-TE1_MEGs_log2FC$Informative_read_log2FC),]
TE1_MEGs_EE_log2FC <- TE1_MEGs_in_EE[,c(1:8,10)]
names(TE1_MEGs_log2FC)[names(TE1_MEGs_log2FC) == "Informative_read_log2FC"] <- "log2FC_TE1"
names(TE1_MEGs_EE_log2FC)[names(TE1_MEGs_EE_log2FC) == "Informative_read_log2FC"] <- "log2FC_EE"
names(TE1_MEGs_log2FC)[names(TE1_MEGs_log2FC) == "Adj.Pvalue"] <- "Adj.P_TE1"
names(TE1_MEGs_EE_log2FC)[names(TE1_MEGs_EE_log2FC) == "Adj.Pvalue"] <- "Adj.P_EE"
TE1_MEGs_in_EE <- full_join(TE1_MEGs_log2FC, TE1_MEGs_EE_log2FC, by="Gene", copy=FALSE,suffix=c("1","2")) 
TE1_MEGs_in_EE$log2FC_EE <- as.character(TE1_MEGs_in_EE$log2FC_EE)
TE1_MEGs_in_EE$log2FC_EE[is.na(TE1_MEGs_in_EE$log2FC_EE)] <- "Not present"

TE1_PEGs_log2FC <- TE1_PEG[,c(1:8,10)]
TE1_PEGs_log2FC <- TE1_PEGs_log2FC[order(-TE1_PEGs_log2FC$Informative_read_log2FC),]
TE1_PEGs_EE_log2FC <- TE1_PEGs_in_EE[,c(1:8,10)]
names(TE1_PEGs_log2FC)[names(TE1_PEGs_log2FC) == "Informative_read_log2FC"] <- "log2FC_TE1"
names(TE1_PEGs_EE_log2FC)[names(TE1_PEGs_EE_log2FC) == "Informative_read_log2FC"] <- "log2FC_EE"
names(TE1_PEGs_log2FC)[names(TE1_PEGs_log2FC) == "Adj.Pvalue"] <- "Adj.P_TE1"
names(TE1_PEGs_EE_log2FC)[names(TE1_PEGs_EE_log2FC) == "Adj.Pvalue"] <- "Adj.P_EE"
TE1_PEGs_in_EE <- full_join(TE1_PEGs_log2FC, TE1_PEGs_EE_log2FC, by="Gene", copy=FALSE,suffix=c("1","2")) 
TE1_PEGs_in_EE$log2FC_EE <- as.character(TE1_PEGs_in_EE$log2FC_EE)
TE1_PEGs_in_EE$log2FC_EE[is.na(TE1_PEGs_in_EE$log2FC_EE)] <- "Not present"

TE1_EE <- rbind(TE1_MEGs_in_EE, TE1_PEGs_in_EE)

#TE1 specific MEGs and PEGs in ESR
TE1_MEGs_log2FC <- TE1_MEG[,c(1:8,10)]
TE1_MEGs_log2FC <- TE1_MEGs_log2FC[order(-TE1_MEGs_log2FC$Informative_read_log2FC),]
TE1_MEGs_ESR_log2FC <- TE1_MEGs_in_ESR[,c(1:8,10)]
names(TE1_MEGs_log2FC)[names(TE1_MEGs_log2FC) == "Informative_read_log2FC"] <- "log2FC_TE1"
names(TE1_MEGs_ESR_log2FC)[names(TE1_MEGs_ESR_log2FC) == "Informative_read_log2FC"] <- "log2FC_ESR"
names(TE1_MEGs_log2FC)[names(TE1_MEGs_log2FC) == "Adj.Pvalue"] <- "Adj.P_TE1"
names(TE1_MEGs_ESR_log2FC)[names(TE1_MEGs_ESR_log2FC) == "Adj.Pvalue"] <- "Adj.P_ESR"
TE1_MEGs_in_ESR <- full_join(TE1_MEGs_log2FC, TE1_MEGs_ESR_log2FC, by="Gene", copy=FALSE,suffix=c("1","2")) 
TE1_MEGs_in_ESR$log2FC_ESR <- as.character(TE1_MEGs_in_ESR$log2FC_ESR)
TE1_MEGs_in_ESR$log2FC_ESR[is.na(TE1_MEGs_in_ESR$log2FC_ESR)] <- "Not present"

TE1_PEGs_log2FC <- TE1_PEG[,c(1:8,10)]
TE1_PEGs_log2FC <- TE1_PEGs_log2FC[order(-TE1_PEGs_log2FC$Informative_read_log2FC),]
TE1_PEGs_ESR_log2FC <- TE1_PEGs_in_ESR[,c(1:8,10)]
names(TE1_PEGs_log2FC)[names(TE1_PEGs_log2FC) == "Informative_read_log2FC"] <- "log2FC_TE1"
names(TE1_PEGs_ESR_log2FC)[names(TE1_PEGs_ESR_log2FC) == "Informative_read_log2FC"] <- "log2FC_ESR"
names(TE1_PEGs_log2FC)[names(TE1_PEGs_log2FC) == "Adj.Pvalue"] <- "Adj.P_TE1"
names(TE1_PEGs_ESR_log2FC)[names(TE1_PEGs_ESR_log2FC) == "Adj.Pvalue"] <- "Adj.P_ESR"
TE1_PEGs_in_ESR <- full_join(TE1_PEGs_log2FC, TE1_PEGs_ESR_log2FC, by="Gene", copy=FALSE,suffix=c("1","2")) 
TE1_PEGs_in_ESR$log2FC_ESR <- as.character(TE1_PEGs_in_ESR$log2FC_ESR)
TE1_PEGs_in_ESR$log2FC_ESR[is.na(TE1_PEGs_in_ESR$log2FC_ESR)] <- "Not present"

TE1_ESR <- rbind(TE1_MEGs_in_ESR, TE1_PEGs_in_ESR)

#Plot ESR MEGs and PEGs as log2FC in TE1 vs log2FC in ESR to find MEGs/PEGs that are different between domains
ESR_all$log2FC_TE1 <- as.numeric(as.character(ESR_all$log2FC_TE1))
ESR_all$log2FC_ESR <- as.numeric(as.character(ESR_all$log2FC_ESR))
TE1_ESR$log2FC_TE1 <- as.numeric(as.character(TE1_ESR$log2FC_TE1))
TE1_ESR$log2FC_ESR <- as.numeric(as.character(TE1_ESR$log2FC_ESR))

#Add a column for the plot for domain specific MEGs/PEGs and overlapping genes to specify which group they belong (TE1, ESR, overlap)
TE1_all_and_ESR_all <- subset(ESR_all, (Gene %in% TE1_ESR$Gene))
TE1_all_and_ESR_all$Domain <- c("ESR+TE1")
TE1_all_only_TE1_ESR <- subset(TE1_ESR, !(Gene %in% ESR_all$Gene))
TE1_all_only_TE1_ESR <- TE1_all_only_TE1_ESR[,c(1,10:17,2:9)]
TE1_all_only_TE1_ESR$Domain <- c("TE1")
ESR_all_only_ESR <- subset(ESR_all, !(Gene %in% TE1_ESR$Gene))
ESR_all_only_ESR$Domain <- c("ESR")

ESR_TE1_all <- rbind(TE1_all_and_ESR_all, TE1_all_only_TE1_EE, ESR_all_only_ESR)

#Plot EE MEGs and PEGs as log2FC in TE1 vs log2FC in EE to find MEGs/PEGs that are different between domains
EE_all$log2FC_TE1 <- as.numeric(as.character(EE_all$log2FC_TE1))
EE_all$log2FC_EE <- as.numeric(as.character(EE_all$log2FC_EE))
TE1_EE$log2FC_TE1 <- as.numeric(as.character(TE1_EE$log2FC_TE1))
TE1_EE$log2FC_EE <- as.numeric(as.character(TE1_EE$log2FC_EE))

#Add a column for the plot for domain specific MEGs/PEGs and overlapping genes to specify which group they belong (TE1, EE, overlap)
TE1_all_and_EE_all <- subset(EE_all, (Gene %in% TE1_EE$Gene))
TE1_all_and_EE_all$Domain <- c("EE+TE1")
TE1_all_only_TE1 <- subset(TE1_EE, !(Gene %in% EE_all$Gene))
TE1_all_only_TE1 <- TE1_all_only_TE1[,c(1,10:17,2:9)]
TE1_all_only_TE1$Domain <- c("TE1")
EE_all_only_EE <- subset(EE_all, !(Gene %in% TE1_EE$Gene))
EE_all_only_EE$Domain <- c("EE")

EE_TE1_all <- rbind(TE1_all_and_EE_all, TE1_all_only_TE1, EE_all_only_EE)

#Filter all samples for >30 informative reads (This should not be necessary, but is extra)
EE_TE1_all <- cbind(EE_TE1_all, All_inf_reads_EE=rowSums(EE_TE1_all[,2:7])) %>% 
  cbind(., All_inf_reads_TE1=rowSums(.[,10:15])) %>%  .[which(.$All_inf_reads_EE > 30),] %>% 
  .[which(.$All_inf_reads_TE1 > 30),] %>% 
  subset(., select=-c(All_inf_reads_EE,All_inf_reads_TE1))

ESR_TE1_all <- cbind(ESR_TE1_all, All_inf_reads_ESR=rowSums(ESR_TE1_all[,2:7])) %>% 
  cbind(., All_inf_reads_TE1=rowSums(.[,10:15])) %>%  .[which(.$All_inf_reads_ESR > 30),] %>% 
  .[which(.$All_inf_reads_TE1 > 30),] %>% 
  subset(., select=-c(All_inf_reads_ESR,All_inf_reads_TE1))

plotp_TE1_EE <- ggplot(EE_TE1_all, aes(x=log2FC_EE, y=log2FC_TE1, shape=Domain)) + 
  geom_point(size=1.5) + 
  scale_shape_manual(values=c(3, 16, 17))+
  xlim(-10, 10) +
  ylim(-10, 10) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") + 
  geom_abline(intercept = -2, slope = 1, linetype="dashed", colour=('#E41A1C')) + 
  geom_abline(intercept = 2, slope = 1, linetype="dashed", colour=('#E41A1C')) +  
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = -0.5) +
  geom_hline(yintercept = 1.5) +
  geom_hline(yintercept = -1.5) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.5) +
  geom_vline(xintercept = -0.5) +
  geom_vline(xintercept = 1.5) +
  geom_vline(xintercept = -1.5) +
  annotate("text", x=-9, y=5, label="TE1: MEG", size=5) +
  annotate("text", x=-9, y=4.5, label="EE: PEG", size=5) +  
  annotate("text", x=9, y=5, label="TE1: MEG", size=5) +
  annotate("text", x=9, y=4.5, label="EE: MEG", size=5) +  
  annotate("text", x=-9, y=-5, label="TE1: PEG", size=5) +
  annotate("text", x=-9, y=-4.5, label="EE: PEG", size=5) +  
  annotate("text", x=9, y=-5, label="TE1: PEG", size=5) +
  annotate("text", x=9, y=-4.5, label="EE: MEG", size=5) +
  theme_classic() + 
  theme(legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(legend.position = c(0.2, 0.9)) + 
  theme(legend.title = element_text(size=10))+ 
  theme(legend.text= element_text(size=10))
plotp_TE1_EE
ggsave(filename= "Figure 5 log2FC_of_TE1_vs_EE.pdf", plot=plotp_TE1_EE, path='./', scale = 1, height = 11.00, width = 11.00)


plotp_TE1_ESR <- ggplot(ESR_TE1_all, aes(x=log2FC_ESR, y=log2FC_TE1, shape=Domain)) + 
  geom_point(size=1.5) + 
  scale_shape_manual(values=c(3, 16, 17)) +
  xlim(-10, 10) +
  ylim(-10, 10) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") + 
  geom_abline(intercept = -2, slope = 1, linetype="dashed", colour=('#E41A1C')) + 
  geom_abline(intercept = 2, slope = 1, linetype="dashed", colour=('#E41A1C')) +  
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = -0.5) +
  geom_hline(yintercept = 1.5) +
  geom_hline(yintercept = -1.5) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.5) +
  geom_vline(xintercept = -0.5) +
  geom_vline(xintercept = 1.5) +
  geom_vline(xintercept = -1.5) +
  annotate("text", x=-9, y=5, label="TE1: MEG", size=5) +
  annotate("text", x=-9, y=4.5, label="ESR: PEG", size=5) +  
  annotate("text", x=9, y=5, label="TE1: MEG", size=5) +
  annotate("text", x=9, y=4.5, label="ESR: MEG", size=5) +  
  annotate("text", x=-9, y=-5, label="TE1: PEG", size=5) +
  annotate("text", x=-9, y=-4.5, label="ESR: PEG", size=5) +  
  annotate("text", x=9, y=-5, label="TE1: PEG", size=5) +
  annotate("text", x=9, y=-4.5, label="ESR: MEG", size=5) +
  theme_classic() + 
  theme(legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(legend.position = c(0.2, 0.9)) + 
  theme(legend.title = element_text(size=10))+ 
  theme(legend.text= element_text(size=10))
plotp_TE1_ESR
ggsave(filename= "Figure 6 log2FC_of_TE1_vs_ESR.pdf", plot=plotp_TE1_ESR, path='.', scale = 1, height = 11.00, width = 11.00)

# List the genes the are removed and not in the plot
EE_TE1_removed <- rbind(TE1_all_and_EE_all, TE1_all_only_TE1, EE_all_only_EE)
ESR_TE1_removed <- rbind(TE1_all_and_ESR_all, TE1_all_only_TE1_ESR, ESR_all_only_ESR)

#Which genes are removed
EE_TE1_removed2 <- cbind(EE_TE1_removed, All_inf_reads_EE=rowSums(EE_TE1_removed[,2:7])) %>% 
  cbind(., All_inf_reads_TE1=rowSums(.[,10:15])) %>% .[rowSums(is.na(.)) > 0, ]  

ESR_TE1_removed2 <- cbind(ESR_TE1_removed, All_inf_reads_ESR=rowSums(ESR_TE1_removed[,2:7])) %>% 
  cbind(., All_inf_reads_TE1=rowSums(.[,10:15])) %>% .[rowSums(is.na(.)) > 0, ] 

write_xlsx(list("EE_and_TE1_MEGs_and_PEGs"=EE_TE1_all, "Removed"=EE_TE1_removed2), "./Output/SData 11 Seed stage specific imprinting.xlsx")
write_xlsx(list("ESR_and_TE1_MEGs_and_PEGs"=ESR_TE1_all, "Removed"=ESR_TE1_removed2), "./Output/SData 12 Domain specific imprinting.xlsx")
