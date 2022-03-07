library(tidyverse)
library(ggplot2)
library(VennDiagram)
library(writexl)
library(readxl)
#R version 4.1.1 (2021-08-10)
#RStudio version 2021.09.0+351 
#tidyverse version 1.3.1
#ggplot2 version 3.3.5
#VennDiagram version 1.6.20
#Writexl version 1.4.0

#Read all tables
EE <- read.csv("./Total inf read and ecotype specific gene filter/EE_inf_reads_homozygous_reads_and_gene_description.csv", header=TRUE)
ESR <- read.csv("./Total inf read and ecotype specific gene filter/ESR_inf_reads_homozygous_reads_and_gene_description.csv", header=TRUE)
TE1 <- read.csv("./Total inf read and ecotype specific gene filter/TE1_inf_reads_homozygous_reads_and_gene_description.csv", header=TRUE)

#MEGs
#Adj.Pvalue < 0.05 and log2FC > 0
EE_MEG <- EE[EE$Informative_read_log2FC > 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
ESR_MEG <- ESR[ESR$Informative_read_log2FC > 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
TE1_MEG <- TE1[TE1$Informative_read_log2FC > 0, ] %>% .[which(.$Adj.Pvalue < 0.05),]

#PEGs
#Adj.Pvalue < 0.05 and log2FC < 0
EE_PEG <- EE[EE$Informative_read_log2FC < 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
ESR_PEG <- ESR[ESR$Informative_read_log2FC < 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
TE1_PEG <- TE1[TE1$Informative_read_log2FC < 0, ] %>% .[which(.$Adj.Pvalue < 0.05),]

#All domains megs and PEGs
Ekelenburg_EE_MEG_gene <- EE_MEG[order(EE_MEG$Gene),] %>% subset(., select=Gene)
Ekelenburg_EE_PEG_gene <- EE_PEG[order(EE_PEG$Gene),] %>% subset(., select=Gene)
Ekelenburg_ESR_MEG_gene <- ESR_MEG[order(ESR_MEG$Gene),] %>% subset(., select=Gene)
Ekelenburg_ESR_PEG_gene <- ESR_PEG[order(ESR_PEG$Gene),] %>% subset(., select=Gene)
Ekelenburg_TE1_MEG_gene <- TE1_MEG[order(TE1_MEG$Gene),] %>% subset(., select=Gene)
Ekelenburg_TE1_PEG_gene <- TE1_PEG[order(TE1_PEG$Gene),] %>% subset(., select=Gene)

All_domains_MEG <- full_join(Ekelenburg_EE_MEG_gene, Ekelenburg_ESR_MEG_gene, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Ekelenburg_TE1_MEG_gene, by="Gene", copy=FALSE,suffix=c("1","2")) 
All_domains_PEG <- full_join(Ekelenburg_EE_PEG_gene, Ekelenburg_ESR_PEG_gene, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Ekelenburg_TE1_PEG_gene, by="Gene", copy=FALSE,suffix=c("1","2"))

#Gene descriptions
Description <- read_excel("./Input/All_gene_descriptions.xlsx")
All_domains_MEG <- subset(Description, Gene %in% All_domains_MEG$Gene) %>% .[order(.$Gene),]
All_domains_PEG <- subset(Description, Gene %in% All_domains_PEG$Gene) %>% .[order(.$Gene),]

#Write MEGs and PEGs for each domain in one excel file 
write_xlsx(list("EE MEG"=EE_MEG, "EE PEG"=EE_PEG, "ESR MEG"=ESR_MEG, "ESR PEG"=ESR_PEG, "TE1 MEG"=TE1_MEG, "TE1 PEG"=TE1_PEG, "All domains MEG"=All_domains_MEG, "All domains PEG"=All_domains_PEG), "./Output/SData 6 All MEGs and PEGs for EE, ESR and TE1.xlsx") 

#To look for overlap, only gene names are necessary
EE_MEG <- EE_MEG[,1, drop=FALSE]
ESR_MEG <- ESR_MEG[,1, drop=FALSE]
TE1_MEG <- TE1_MEG[,1, drop=FALSE]
EE_PEG <- EE_PEG[,1, drop=FALSE]
ESR_PEG <- ESR_PEG[,1, drop=FALSE]
TE1_PEG <- TE1_PEG[,1, drop=FALSE]

#Domain specific MEGs and PEGs
EE_MEG_specific <- subset(EE_MEG, !(Gene %in% ESR_MEG$Gene) & !(Gene %in% TE1_MEG$Gene)) 
ESR_MEG_specific <- subset(ESR_MEG, !(Gene %in% EE_MEG$Gene) & !(Gene %in% TE1_MEG$Gene))
TE1_MEG_specific <- subset(TE1_MEG, !(Gene %in% EE_MEG$Gene) & !(Gene %in% ESR_MEG$Gene))

EE_PEG_specific <- subset(EE_PEG, !(Gene %in% ESR_PEG$Gene) & !(Gene %in% TE1_PEG$Gene)) 
ESR_PEG_specific <- subset(ESR_PEG, !(Gene %in% EE_PEG$Gene) & !(Gene %in% TE1_PEG$Gene))
TE1_PEG_specific <- subset(TE1_PEG, !(Gene %in% EE_PEG$Gene) & !(Gene %in% ESR_PEG$Gene))

#Overlapping MEGs and PEGs
EE_ESR_MEG <- subset(EE_MEG, Gene %in% ESR_MEG$Gene & !(Gene %in% TE1_MEG$Gene))
EE_TE1_MEG <- subset(EE_MEG, Gene %in% TE1_MEG$Gene & !(Gene %in% ESR_MEG$Gene))
ESR_TE1_MEG <- subset(ESR_MEG, Gene %in% TE1_MEG$Gene & !(Gene %in% EE_MEG$Gene))
EE_ESR_TE1_MEG <- subset(EE_MEG, Gene %in% ESR_MEG$Gene) %>% subset(., Gene %in% TE1_MEG$Gene) 

EE_ESR_PEG <- subset(EE_PEG, Gene %in% ESR_PEG$Gene & !(Gene %in% TE1_PEG$Gene))
EE_TE1_PEG <- subset(EE_PEG, Gene %in% TE1_PEG$Gene & !(Gene %in% ESR_PEG$Gene))
ESR_TE1_PEG <- subset(ESR_PEG, Gene %in% TE1_PEG$Gene & !(Gene %in% EE_PEG$Gene))
EE_ESR_TE1_PEG <- subset(EE_PEG, Gene %in% ESR_PEG$Gene) %>% subset(., Gene %in% TE1_PEG$Gene) 

#Write MEGs and PEGs for each domain in one excel file 
write_xlsx(list("EE specific MEGs"=EE_MEG_specific, "EE specific PEGs"=EE_PEG_specific, "ESR specific MEGs"=ESR_MEG_specific, "ESR specific PEGs"=ESR_PEG_specific, "TE1 specific MEGs"=TE1_MEG_specific, "TE1 specific PEGs"=TE1_PEG_specific, "EE and ESR MEGs"=EE_ESR_MEG, "EE and ESR PEGs"=EE_ESR_PEG, "EE and TE1 MEGs"=EE_TE1_MEG, "EE and TE1 PEGs"=EE_TE1_PEG, "ESR and TE1 MEGs"=ESR_TE1_MEG, "ESR and TE1 PEGs"=ESR_TE1_PEG, "EE, ESR and TE1 MEGs"=EE_ESR_TE1_MEG, "EE, ESR and TE1 PEGs"=EE_ESR_TE1_PEG), "./Output/SData 10 Domains specific and overlapping MEGs and PEGs.xlsx")
