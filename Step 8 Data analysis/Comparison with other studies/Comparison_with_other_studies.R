library(tidyverse)
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(writexl)
library(readxl)

#R version 4.1.1 (2021-08-10)
#RStudio version 2021.09.0+351 
#tidyverse version 1.3.1 
#ggplot2 version 3.3.5
#VennDiagram version 1.6.20
#Writexl version 1.4.0
#dplyr version 1.0.7

#Read all tables
EE <- read.csv("./EEpos_no_transcript.csv", header=TRUE)
ESR <- read.csv("./ESRpos_no_transcript.csv", header=TRUE)
TE1 <- read.csv("./TE1pos_no_transcript.csv", header=TRUE)

Other <- read_tsv("./Other_studies2.tsv") %>% na_if("") 

#List of genes that are filtered for ecotype specific filter
Total <- full_join(EE, ESR, by="Gene", copy=FALSE,suffix=c("1","2"))  %>% full_join(., TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 

#Ekelenburg et al MEGs filtered
Ekelenburg_EE_MEG <- EE[EE$Informative_read_log2FC >= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
Ekelenburg_ESR_MEG <- ESR[ESR$Informative_read_log2FC >= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
Ekelenburg_TE1_MEG <- TE1[TE1$Informative_read_log2FC >= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 

#Ekelenburg et al PEGs filtered
Ekelenburg_EE_PEG <- EE[EE$Informative_read_log2FC <= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
Ekelenburg_ESR_PEG <- ESR[ESR$Informative_read_log2FC <= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
Ekelenburg_TE1_PEG <- TE1[TE1$Informative_read_log2FC <= 0, ] %>% .[which(.$Adj.Pvalue < 0.05),] 

#Generate lists of MEGs and PEGs of all domains together
Ekelenburg_EE_MEG_gene <- Ekelenburg_EE_MEG[order(Ekelenburg_EE_MEG$Gene),] %>% subset(., select=Gene)
Ekelenburg_EE_PEG_gene <- Ekelenburg_EE_PEG[order(Ekelenburg_EE_PEG$Gene),] %>% subset(., select=Gene)
Ekelenburg_ESR_MEG_gene <- Ekelenburg_ESR_MEG[order(Ekelenburg_ESR_MEG$Gene),] %>% subset(., select=Gene)
Ekelenburg_ESR_PEG_gene <- Ekelenburg_ESR_PEG[order(Ekelenburg_ESR_PEG$Gene),] %>% subset(., select=Gene)
Ekelenburg_TE1_MEG_gene <- Ekelenburg_TE1_MEG[order(Ekelenburg_TE1_MEG$Gene),] %>% subset(., select=Gene)
Ekelenburg_TE1_PEG_gene <- Ekelenburg_TE1_PEG[order(Ekelenburg_TE1_PEG$Gene),] %>% subset(., select=Gene)

MEGs <- full_join(Ekelenburg_EE_MEG_gene, Ekelenburg_ESR_MEG_gene, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Ekelenburg_TE1_MEG_gene, by="Gene", copy=FALSE,suffix=c("1","2")) 
PEGs <- full_join(Ekelenburg_EE_PEG_gene, Ekelenburg_ESR_PEG_gene, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Ekelenburg_TE1_PEG_gene, by="Gene", copy=FALSE,suffix=c("1","2")) 

#Export lists
write_xlsx(list("EE_MEG"=Ekelenburg_EE_MEG_gene, "EE_PEG"=Ekelenburg_EE_PEG_gene, "ESR_MEG"=Ekelenburg_ESR_MEG_gene, "ESR_PEG"=Ekelenburg_ESR_PEG_gene, "TE1_MEG"=Ekelenburg_TE1_MEG_gene, "TE1_PEG"=Ekelenburg_TE1_PEG_gene, "All_MEG"=MEGs, "All_PEG"=PEGs), "./All_domains_MEGs_and_PEGs_without_transcriptnumber.xlsx")

#Overlap between Pignatta et al and van Ekelenburg et al
#EE
Pignatta_MEGs <- subset(Other, select=c(Pignatta.et.al.MEGs)) %>% drop_na()
Pignatta_MEGs_in_EE <- subset(EE, (Gene %in% Pignatta_MEGs$Pignatta.et.al.MEGs))
Pignatta_PEGs <- subset(Other, select=c(Pignatta.et.al.PEGs)) %>% drop_na()
Pignatta_PEGs_in_EE  <- subset(EE, (Gene %in% Pignatta_PEGs$Pignatta.et.al.PEGs))

Ekelenburg_Pignatta_EE_MEGs  <- subset(Ekelenburg_EE_MEG, (Gene %in% Pignatta_MEGs$Pignatta.et.al.MEGs))
Ekelenburg_Pignatta_EE_PEGs  <- subset(Ekelenburg_EE_PEG, (Gene %in% Pignatta_PEGs$Pignatta.et.al.PEGs))

#ESR
Pignatta_MEGs <- subset(Other, select=c(Pignatta.et.al.MEGs)) %>% drop_na()
Pignatta_MEGs_in_ESR <- subset(ESR, (Gene %in% Pignatta_MEGs$Pignatta.et.al.MEGs))
Pignatta_PEGs <- subset(Other, select=c(Pignatta.et.al.PEGs)) %>% drop_na()
Pignatta_PEGs_in_ESR  <- subset(ESR, (Gene %in% Pignatta_PEGs$Pignatta.et.al.PEGs))

Ekelenburg_Pignatta_ESR_MEGs  <- subset(Ekelenburg_ESR_MEG, (Gene %in% Pignatta_MEGs$Pignatta.et.al.MEGs))
Ekelenburg_Pignatta_ESR_PEGs  <- subset(Ekelenburg_ESR_PEG, (Gene %in% Pignatta_PEGs$Pignatta.et.al.PEGs))

#TE1
Pignatta_MEGs <- subset(Other, select=c(Pignatta.et.al.MEGs)) %>% drop_na()
Pignatta_MEGs_in_TE1 <- subset(TE1, (Gene %in% Pignatta_MEGs$Pignatta.et.al.MEGs))
Pignatta_PEGs <- subset(Other, select=c(Pignatta.et.al.PEGs)) %>% drop_na()
Pignatta_PEGs_in_TE1  <- subset(TE1, (Gene %in% Pignatta_PEGs$Pignatta.et.al.PEGs))

Ekelenburg_Pignatta_TE1_MEGs  <- subset(Ekelenburg_TE1_MEG, (Gene %in% Pignatta_MEGs$Pignatta.et.al.MEGs))
Ekelenburg_Pignatta_TE1_PEGs  <- subset(Ekelenburg_TE1_PEG, (Gene %in% Pignatta_PEGs$Pignatta.et.al.PEGs))

#Overlap between Picard et al and van Ekelenburg et al
#EE
Picard_MEGs <- subset(Other, select=c(Picard.et.al.MEGs)) %>% drop_na()
Picard_MEGs_in_EE <- subset(EE, (Gene %in% Picard_MEGs$Picard.et.al.MEGs))
Picard_PEGs <- subset(Other, select=c(Picard.et.al.PEGs)) %>% drop_na()
Picard_PEGs_in_EE  <- subset(EE, (Gene %in% Picard_PEGs$Picard.et.al.PEGs))

Ekelenburg_Picard_EE_MEGs  <- subset(Ekelenburg_EE_MEG, (Gene %in% Picard_MEGs$Picard.et.al.MEGs))
Ekelenburg_Picard_EE_PEGs  <- subset(Ekelenburg_EE_PEG, (Gene %in% Picard_PEGs$Picard.et.al.PEGs))

#ESR
Picard_MEGs <- subset(Other, select=c(Picard.et.al.MEGs)) %>% drop_na()
Picard_MEGs_in_ESR <- subset(ESR, (Gene %in% Picard_MEGs$Picard.et.al.MEGs))
Picard_PEGs <- subset(Other, select=c(Picard.et.al.PEGs)) %>% drop_na()
Picard_PEGs_in_ESR  <- subset(ESR, (Gene %in% Picard_PEGs$Picard.et.al.PEGs))

Ekelenburg_Picard_ESR_MEGs  <- subset(Ekelenburg_ESR_MEG, (Gene %in% Picard_MEGs$Picard.et.al.MEGs))
Ekelenburg_Picard_ESR_PEGs  <- subset(Ekelenburg_ESR_PEG, (Gene %in% Picard_PEGs$Picard.et.al.PEGs))

#TE1
Picard_MEGs <- subset(Other, select=c(Picard.et.al.MEGs)) %>% drop_na()
Picard_MEGs_in_TE1 <- subset(TE1, (Gene %in% Picard_MEGs$Picard.et.al.MEGs))
Picard_PEGs <- subset(Other, select=c(Picard.et.al.PEGs)) %>% drop_na()
Picard_PEGs_in_TE1  <- subset(TE1, (Gene %in% Picard_PEGs$Picard.et.al.PEGs))

Ekelenburg_Picard_TE1_MEGs  <- subset(Ekelenburg_TE1_MEG, (Gene %in% Picard_MEGs$Picard.et.al.MEGs))
Ekelenburg_Picard_TE1_PEGs  <- subset(Ekelenburg_TE1_PEG, (Gene %in% Picard_PEGs$Picard.et.al.PEGs))

#Overlap between Hornslien et al and van Ekelenburg et al
#EE
Hornslien_MEGs <- subset(Other, select=c(Hornslien.et.al.MEGs)) %>% drop_na()
Hornslien_MEGs_in_EE <- subset(EE, (Gene %in% Hornslien_MEGs$Hornslien.et.al.MEGs))
Hornslien_PEGs <- subset(Other, select=c(Hornslien.et.al.PEGs)) %>% drop_na()
Hornslien_PEGs_in_EE  <- subset(EE, (Gene %in% Hornslien_PEGs$Hornslien.et.al.PEGs))

Ekelenburg_Hornslien_EE_MEGs  <- subset(Ekelenburg_EE_MEG, (Gene %in% Hornslien_MEGs$Hornslien.et.al.MEGs))
Ekelenburg_Hornslien_EE_PEGs  <- subset(Ekelenburg_EE_PEG, (Gene %in% Hornslien_PEGs$Hornslien.et.al.PEGs))

#ESR
Hornslien_MEGs <- subset(Other, select=c(Hornslien.et.al.MEGs)) %>% drop_na()
Hornslien_MEGs_in_ESR <- subset(ESR, (Gene %in% Hornslien_MEGs$Hornslien.et.al.MEGs))
Hornslien_PEGs <- subset(Other, select=c(Hornslien.et.al.PEGs)) %>% drop_na()
Hornslien_PEGs_in_ESR  <- subset(ESR, (Gene %in% Hornslien_PEGs$Hornslien.et.al.PEGs))

Ekelenburg_Hornslien_ESR_MEGs  <- subset(Ekelenburg_ESR_MEG, (Gene %in% Hornslien_MEGs$Hornslien.et.al.MEGs))
Ekelenburg_Hornslien_ESR_PEGs  <- subset(Ekelenburg_ESR_PEG, (Gene %in% Hornslien_PEGs$Hornslien.et.al.PEGs))

#TE1
Hornslien_MEGs <- subset(Other, select=c(Hornslien.et.al.MEGs)) %>% drop_na()
Hornslien_MEGs_in_TE1 <- subset(TE1, (Gene %in% Hornslien_MEGs$Hornslien.et.al.MEGs))
Hornslien_PEGs <- subset(Other, select=c(Hornslien.et.al.PEGs)) %>% drop_na()
Hornslien_PEGs_in_TE1  <- subset(TE1, (Gene %in% Hornslien_PEGs$Hornslien.et.al.PEGs))

Ekelenburg_Hornslien_TE1_MEGs  <- subset(Ekelenburg_TE1_MEG, (Gene %in% Hornslien_MEGs$Hornslien.et.al.MEGs))
Ekelenburg_Hornslien_TE1_PEGs  <- subset(Ekelenburg_TE1_PEG, (Gene %in% Hornslien_PEGs$Hornslien.et.al.PEGs))

#Overlap between Del_Toro et al and van Ekelenburg et al
#EE
Del_Toro_MEGs <- subset(Other, select=c(Del_Toro.et.al.MEGs)) %>% drop_na()
Del_Toro_MEGs_in_EE <- subset(EE, (Gene %in% Del_Toro_MEGs$Del_Toro.et.al.MEGs))
Del_Toro_PEGs <- subset(Other, select=c(Del_Toro.et.al.PEGs)) %>% drop_na()
Del_Toro_PEGs_in_EE  <- subset(EE, (Gene %in% Del_Toro_PEGs$Del_Toro.et.al.PEGs))

Ekelenburg_Del_Toro_EE_MEGs  <- subset(Ekelenburg_EE_MEG, (Gene %in% Del_Toro_MEGs$Del_Toro.et.al.MEGs))
Ekelenburg_Del_Toro_EE_PEGs  <- subset(Ekelenburg_EE_PEG, (Gene %in% Del_Toro_PEGs$Del_Toro.et.al.PEGs))

#ESR
Del_Toro_MEGs <- subset(Other, select=c(Del_Toro.et.al.MEGs)) %>% drop_na()
Del_Toro_MEGs_in_ESR <- subset(ESR, (Gene %in% Del_Toro_MEGs$Del_Toro.et.al.MEGs))
Del_Toro_PEGs <- subset(Other, select=c(Del_Toro.et.al.PEGs)) %>% drop_na()
Del_Toro_PEGs_in_ESR  <- subset(ESR, (Gene %in% Del_Toro_PEGs$Del_Toro.et.al.PEGs))

Ekelenburg_Del_Toro_ESR_MEGs  <- subset(Ekelenburg_ESR_MEG, (Gene %in% Del_Toro_MEGs$Del_Toro.et.al.MEGs))
Ekelenburg_Del_Toro_ESR_PEGs  <- subset(Ekelenburg_ESR_PEG, (Gene %in% Del_Toro_PEGs$Del_Toro.et.al.PEGs))

#TE1
Del_Toro_MEGs <- subset(Other, select=c(Del_Toro.et.al.MEGs)) %>% drop_na()
Del_Toro_MEGs_in_TE1 <- subset(TE1, (Gene %in% Del_Toro_MEGs$Del_Toro.et.al.MEGs))
Del_Toro_PEGs <- subset(Other, select=c(Del_Toro.et.al.PEGs)) %>% drop_na()
Del_Toro_PEGs_in_TE1  <- subset(TE1, (Gene %in% Del_Toro_PEGs$Del_Toro.et.al.PEGs))

Ekelenburg_Del_Toro_TE1_MEGs  <- subset(Ekelenburg_TE1_MEG, (Gene %in% Del_Toro_MEGs$Del_Toro.et.al.MEGs))
Ekelenburg_Del_Toro_TE1_PEGs  <- subset(Ekelenburg_TE1_PEG, (Gene %in% Del_Toro_PEGs$Del_Toro.et.al.PEGs))

#Genes of the other studies that are present in all domains together
Hornslien_MEGs_present <- full_join(Hornslien_MEGs_in_ESR, Hornslien_MEGs_in_EE, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Hornslien_MEGs_in_TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 
Hornslien_PEGs_present <- full_join(Hornslien_PEGs_in_ESR, Hornslien_PEGs_in_EE, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Hornslien_PEGs_in_TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 
Pignatta_MEGs_present <- full_join(Pignatta_MEGs_in_ESR, Pignatta_MEGs_in_EE, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Pignatta_MEGs_in_TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 
Pignatta_PEGs_present <- full_join(Pignatta_PEGs_in_ESR, Pignatta_PEGs_in_EE, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Pignatta_PEGs_in_TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 
Picard_MEGs_present <- full_join(Picard_MEGs_in_ESR, Picard_MEGs_in_EE, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Picard_MEGs_in_TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 
Picard_PEGs_present <- full_join(Picard_PEGs_in_ESR, Picard_PEGs_in_EE, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Picard_PEGs_in_TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 
Del_Toro_MEGs_present <- full_join(Del_Toro_MEGs_in_ESR, Del_Toro_MEGs_in_EE, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Del_Toro_MEGs_in_TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 
Del_Toro_PEGs_present <- full_join(Del_Toro_PEGs_in_ESR, Del_Toro_PEGs_in_EE, by="Gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., Del_Toro_PEGs_in_TE1, by="Gene", copy=FALSE,suffix=c("1","2")) 

#Overlap with other studies of All domain MEGs and PEGs 
Pignatta_MEGs_in_All <- subset(MEGs, (Gene %in% Pignatta_MEGs$Pignatta.et.al.MEGs))
Pignatta_PEGs_in_All <- subset(PEGs, (Gene %in% Pignatta_PEGs$Pignatta.et.al.PEGs))
Hornslien_MEGs_in_All <- subset(MEGs, (Gene %in% Hornslien_MEGs$Hornslien.et.al.MEGs))
Hornslien_PEGs_in_All <- subset(PEGs, (Gene %in% Hornslien_PEGs$Hornslien.et.al.PEGs))
Picard_MEGs_in_All <- subset(MEGs, (Gene %in% Picard_MEGs$Picard.et.al.MEGs))
Picard_PEGs_in_All <- subset(PEGs, (Gene %in% Picard_PEGs$Picard.et.al.PEGs))
Del_Toro_MEGs_in_All <- subset(MEGs, (Gene %in% Del_Toro_MEGs$Del_Toro.et.al.MEGs))
Del_Toro_PEGs_in_All <- subset(PEGs, (Gene %in% Del_Toro_PEGs$Del_Toro.et.al.PEGs))

All_domain_MEGs_in_other_studies <- full_join(Pignatta_MEGs_in_All, Hornslien_MEGs_in_All, by="Gene", copy=FALSE,suffix=c("1","2")) %>%
  full_join(., Picard_MEGs_in_All, by="Gene", copy=FALSE,suffix=c("1","2")) %>% 
  full_join(., Del_Toro_MEGs_in_All, by="Gene", copy=FALSE,suffix=c("1","2")) 
All_domain_PEGs_in_other_studies <- full_join(Pignatta_PEGs_in_All, Hornslien_PEGs_in_All, by="Gene", copy=FALSE,suffix=c("1","2")) %>%
  full_join(., Picard_PEGs_in_All, by="Gene", copy=FALSE,suffix=c("1","2")) %>% 
  full_join(., Del_Toro_PEGs_in_All, by="Gene", copy=FALSE,suffix=c("1","2")) 

# Add gene description
Description <- read_excel("./Input/All_gene_descriptions.xlsx")
All_domain_MEGs_in_other_studies <- subset(Description, Gene %in% All_domain_MEGs_in_other_studies$Gene) %>% .[order(.$Gene),]
All_domain_PEGs_in_other_studies <- subset(Description, Gene %in% All_domain_PEGs_in_other_studies$Gene) %>% .[order(.$Gene),]

#Export lists
write_xlsx(list("MEGs"=All_domain_MEGs_in_other_studies, "PEGs"=All_domain_PEGs_in_other_studies), "./SData 8 All domains overlapping MEGs and PEGs with other studies with description.xlsx")

#Compare MEGs and PEGs without filtering other studies
#Renaming columns of each study set
Pignatta_MEGs <- subset(Other, select=c(Pignatta.et.al.MEGs)) %>% drop_na() 
names(Pignatta_MEGs)[names(Pignatta_MEGs) == "Pignatta.et.al.MEGs"] <- "Gene"
Pignatta_PEGs <- subset(Other, select=c(Pignatta.et.al.PEGs)) %>% drop_na() 
names(Pignatta_PEGs)[names(Pignatta_PEGs) == "Pignatta.et.al.PEGs"] <- "Gene"
Hornslien_MEGs <- subset(Other, select=c(Hornslien.et.al.MEGs)) %>% drop_na() 
names(Hornslien_MEGs)[names(Hornslien_MEGs) == "Hornslien.et.al.MEGs"] <- "Gene"
Hornslien_PEGs <- subset(Other, select=c(Hornslien.et.al.PEGs)) %>% drop_na() 
names(Hornslien_PEGs)[names(Hornslien_PEGs) == "Hornslien.et.al.PEGs"] <- "Gene"
Picard_MEGs <- subset(Other, select=c(Picard.et.al.MEGs)) %>% drop_na() 
names(Picard_MEGs)[names(Picard_MEGs) == "Picard.et.al.MEGs"] <- "Gene"
Picard_PEGs <- subset(Other, select=c(Picard.et.al.PEGs)) %>% drop_na() 
names(Picard_PEGs)[names(Picard_PEGs) == "Picard.et.al.PEGs"] <- "Gene"
Del_Toro_MEGs <- subset(Other, select=c(Del_Toro.et.al.MEGs)) %>% drop_na() 
names(Del_Toro_MEGs)[names(Del_Toro_MEGs) == "Del_Toro.et.al.MEGs"] <- "Gene"
Del_Toro_PEGs <- subset(Other, select=c(Del_Toro.et.al.PEGs)) %>% drop_na() 
names(Del_Toro_PEGs)[names(Del_Toro_PEGs) == "Del_Toro.et.al.PEGs"] <- "Gene"

#Table overview showing the number of genes for each dataset and how many genes were filtered out
#MEGs and PEGs for our study
Ekelenburg_ESR_MEGs <- nrow(Ekelenburg_ESR_MEG)
Ekelenburg_ESR_PEGs <- nrow(Ekelenburg_ESR_PEG)
Ekelenburg_EE_MEGs <- nrow(Ekelenburg_EE_MEG)
Ekelenburg_EE_PEGs <- nrow(Ekelenburg_EE_PEG)
Ekelenburg_TE1_MEGs <- nrow(Ekelenburg_TE1_MEG)
Ekelenburg_TE1_PEGs <- nrow(Ekelenburg_TE1_PEG)
Ekelenburg_All_MEGs <- nrow(MEGs)
Ekelenburg_All_PEGs <- nrow(PEGs)

#Number of genes in each study
Picard_All_MEGs <- nrow(Picard_MEGs)
Picard_All_PEGs <- nrow(Picard_PEGs)
Hornslien_All_MEGs <- nrow(Hornslien_MEGs)
Hornslien_All_PEGs <- nrow(Hornslien_PEGs)
Pignatta_All_MEGs <- nrow(Pignatta_MEGs)
Pignatta_All_PEGs <- nrow(Pignatta_PEGs)
Del_Toro_All_MEGs <- nrow(Del_Toro_MEGs)
Del_Toro_All_PEGs <- nrow(Del_Toro_PEGs)

#Number of genes in EE
Hornslien_MEG_present_in_EE <- nrow(Hornslien_MEGs_in_EE)
Hornslien_PEG_present_in_EE <- nrow(Hornslien_PEGs_in_EE)
Picard_MEG_present_in_EE <- nrow(Picard_MEGs_in_EE)
Picard_PEG_present_in_EE <- nrow(Picard_PEGs_in_EE)
Pignatta_MEG_present_in_EE <- nrow(Pignatta_MEGs_in_EE)
Pignatta_PEG_present_in_EE <- nrow(Pignatta_PEGs_in_EE)
Del_Toro_MEG_present_in_EE <- nrow(Del_Toro_MEGs_in_EE)
Del_Toro_PEG_present_in_EE <- nrow(Del_Toro_PEGs_in_EE)

#Number of overlapping genes in EE
Pignatta_EE_MEG_overlap <- nrow(Ekelenburg_Pignatta_EE_MEGs)
Pignatta_EE_PEG_overlap <- nrow(Ekelenburg_Pignatta_EE_PEGs)
Hornslien_EE_MEG_overlap <- nrow(Ekelenburg_Hornslien_EE_MEGs)
Hornslien_EE_PEG_overlap <- nrow(Ekelenburg_Hornslien_EE_PEGs)
Picard_EE_MEG_overlap <- nrow(Ekelenburg_Picard_EE_MEGs)
Picard_EE_PEG_overlap <- nrow(Ekelenburg_Picard_EE_PEGs)
Del_Toro_EE_MEG_overlap <- nrow(Ekelenburg_Del_Toro_EE_MEGs)
Del_Toro_EE_PEG_overlap <- nrow(Ekelenburg_Del_Toro_EE_PEGs)

#Number of genes in ESR
Hornslien_MEG_present_in_ESR <- nrow(Hornslien_MEGs_in_ESR)
Hornslien_PEG_present_in_ESR <- nrow(Hornslien_PEGs_in_ESR)
Picard_MEG_present_in_ESR <- nrow(Picard_MEGs_in_ESR)
Picard_PEG_present_in_ESR <- nrow(Picard_PEGs_in_ESR)
Pignatta_MEG_present_in_ESR <- nrow(Pignatta_MEGs_in_ESR)
Pignatta_PEG_present_in_ESR <- nrow(Pignatta_PEGs_in_ESR)
Del_Toro_MEG_present_in_ESR <- nrow(Del_Toro_MEGs_in_ESR)
Del_Toro_PEG_present_in_ESR <- nrow(Del_Toro_PEGs_in_ESR)

#Number of overlapping genes in ESR
Pignatta_ESR_MEG_overlap <- nrow(Ekelenburg_Pignatta_ESR_MEGs)
Pignatta_ESR_PEG_overlap <- nrow(Ekelenburg_Pignatta_ESR_PEGs)
Hornslien_ESR_MEG_overlap <- nrow(Ekelenburg_Hornslien_ESR_MEGs)
Hornslien_ESR_PEG_overlap <- nrow(Ekelenburg_Hornslien_ESR_PEGs)
Picard_ESR_MEG_overlap <- nrow(Ekelenburg_Picard_ESR_MEGs)
Picard_ESR_PEG_overlap <- nrow(Ekelenburg_Picard_ESR_PEGs)
Del_Toro_ESR_MEG_overlap <- nrow(Ekelenburg_Del_Toro_ESR_MEGs)
Del_Toro_ESR_PEG_overlap <- nrow(Ekelenburg_Del_Toro_ESR_PEGs)

#Number of genes in TE1
Hornslien_MEG_present_in_TE1 <- nrow(Hornslien_MEGs_in_TE1)
Hornslien_PEG_present_in_TE1 <- nrow(Hornslien_PEGs_in_TE1)
Picard_MEG_present_in_TE1 <- nrow(Picard_MEGs_in_TE1)
Picard_PEG_present_in_TE1 <- nrow(Picard_PEGs_in_TE1)
Pignatta_MEG_present_in_TE1 <- nrow(Pignatta_MEGs_in_TE1)
Pignatta_PEG_present_in_TE1 <- nrow(Pignatta_PEGs_in_TE1)
Del_Toro_MEG_present_in_TE1 <- nrow(Del_Toro_MEGs_in_TE1)
Del_Toro_PEG_present_in_TE1 <- nrow(Del_Toro_PEGs_in_TE1)

#Number of overlapping genes in TE1
Pignatta_TE1_MEG_overlap <- nrow(Ekelenburg_Pignatta_TE1_MEGs)
Pignatta_TE1_PEG_overlap <- nrow(Ekelenburg_Pignatta_TE1_PEGs)
Hornslien_TE1_MEG_overlap <- nrow(Ekelenburg_Hornslien_TE1_MEGs)
Hornslien_TE1_PEG_overlap <- nrow(Ekelenburg_Hornslien_TE1_PEGs)
Picard_TE1_MEG_overlap <- nrow(Ekelenburg_Picard_TE1_MEGs)
Picard_TE1_PEG_overlap <- nrow(Ekelenburg_Picard_TE1_PEGs)
Del_Toro_TE1_MEG_overlap <- nrow(Ekelenburg_Del_Toro_TE1_MEGs)
Del_Toro_TE1_PEG_overlap <- nrow(Ekelenburg_Del_Toro_TE1_PEGs)

#Number of genes in all domains
Pignatta_MEG_present_in_All <- nrow(Pignatta_MEGs_present)
Pignatta_PEG_present_in_All <- nrow(Pignatta_PEGs_present)
Hornslien_MEG_present_in_All <- nrow(Hornslien_MEGs_present)
Hornslien_PEG_present_in_All <- nrow(Hornslien_PEGs_present)
Picard_MEG_present_in_All <- nrow(Picard_MEGs_present)
Picard_PEG_present_in_All <- nrow(Picard_PEGs_present)
Del_Toro_MEG_present_in_All <- nrow(Del_Toro_MEGs_present)
Del_Toro_PEG_present_in_All <- nrow(Del_Toro_PEGs_present)

#Number of overlapping genes in All domains
Hornslien_All_MEG_overlap <- nrow(Hornslien_MEGs_in_All)
Hornslien_All_PEG_overlap <- nrow(Hornslien_PEGs_in_All)
Picard_All_MEG_overlap <- nrow(Picard_MEGs_in_All)
Picard_All_PEG_overlap <- nrow(Picard_PEGs_in_All)
Pignatta_All_MEG_overlap <- nrow(Pignatta_MEGs_in_All)
Pignatta_All_PEG_overlap <- nrow(Pignatta_PEGs_in_All)
Del_Toro_All_MEG_overlap <- nrow(Del_Toro_MEGs_in_All)
Del_Toro_All_PEG_overlap <- nrow(Del_Toro_PEGs_in_All)

Study=c('Picard_MEGs', 'Picard_PEGs', 'Del_Toro_MEGs', 'Del_Toro_PEGs', 'Hornslien_MEGs',
        'Hornslien_PEGs', 'Pignatta_MEGs', 'Pignatta_PEGs')
Genes=c(Picard_All_MEGs, Picard_All_PEGs, Del_Toro_All_MEGs, Del_Toro_All_PEGs, Hornslien_All_MEGs, 
        Hornslien_All_PEGs, Pignatta_All_MEGs, Pignatta_All_PEGs)
Filtered_EE=c(Picard_MEG_present_in_EE, Picard_PEG_present_in_EE, Del_Toro_MEG_present_in_EE, Del_Toro_PEG_present_in_EE, 
              Hornslien_MEG_present_in_EE, Hornslien_PEG_present_in_EE, 
              Pignatta_MEG_present_in_EE, Pignatta_PEG_present_in_EE)
Overlap_EE=c(Picard_EE_MEG_overlap, Picard_EE_PEG_overlap, Del_Toro_EE_MEG_overlap, Del_Toro_EE_PEG_overlap, 
             Hornslien_EE_MEG_overlap, Hornslien_EE_PEG_overlap, 
             Pignatta_EE_MEG_overlap, Pignatta_EE_PEG_overlap)
Filtered_ESR=c(Picard_MEG_present_in_ESR, Picard_PEG_present_in_ESR, Del_Toro_MEG_present_in_ESR, Del_Toro_PEG_present_in_ESR, 
               Hornslien_MEG_present_in_ESR, Hornslien_PEG_present_in_ESR, 
               Pignatta_MEG_present_in_ESR, Pignatta_PEG_present_in_ESR)
Overlap_ESR=c(Picard_ESR_MEG_overlap, Picard_ESR_PEG_overlap, Del_Toro_ESR_MEG_overlap, Del_Toro_ESR_PEG_overlap, 
              Hornslien_ESR_MEG_overlap, Hornslien_ESR_PEG_overlap, 
              Pignatta_ESR_MEG_overlap, Pignatta_ESR_PEG_overlap)
Filtered_TE1=c(Picard_MEG_present_in_TE1, Picard_PEG_present_in_TE1, Del_Toro_MEG_present_in_TE1, Del_Toro_PEG_present_in_TE1, 
               Hornslien_MEG_present_in_TE1, Hornslien_PEG_present_in_TE1, 
               Pignatta_MEG_present_in_TE1, Pignatta_PEG_present_in_TE1)
Overlap_TE1=c(Picard_TE1_MEG_overlap, Picard_TE1_PEG_overlap, Del_Toro_TE1_MEG_overlap, Del_Toro_TE1_PEG_overlap, 
              Hornslien_TE1_MEG_overlap, Hornslien_TE1_PEG_overlap,
              Pignatta_TE1_MEG_overlap, Pignatta_TE1_PEG_overlap)
Filtered_All=c(Picard_MEG_present_in_All, Picard_PEG_present_in_All, Del_Toro_MEG_present_in_All, Del_Toro_PEG_present_in_All, 
               Hornslien_MEG_present_in_All, Hornslien_PEG_present_in_All, 
               Pignatta_MEG_present_in_All, Pignatta_PEG_present_in_All)
Overlap_All=c(Picard_All_MEG_overlap, Picard_All_PEG_overlap, Del_Toro_All_MEG_overlap, Del_Toro_All_PEG_overlap, 
              Hornslien_All_MEG_overlap, Hornslien_All_PEG_overlap, 
              Pignatta_All_MEG_overlap, Pignatta_All_PEG_overlap)

Data_EE <- data.frame("Study"=Study, "Identified"=Genes, 
                      "Overlapping"=Overlap_EE, "Overlap percentage"=(Overlap_EE/Genes), 
                      "Analyzed in this study"=Filtered_EE, 
                      "Overlap percentage only analyzed genes"=(Overlap_EE/Filtered_EE))
Data_ESR <- data.frame("Study"=Study, "Identified"=Genes, 
                       "Overlapping"=Overlap_ESR, "Overlap percentage"=(Overlap_ESR/Genes), 
                       "Analyzed in this study"=Filtered_ESR, 
                       "Overlap percentage only analyzed genes"=(Overlap_ESR/Filtered_ESR))
Data_TE1 <- data.frame("Study"=Study, "Identified"=Genes, 
                       "Overlappings"=Overlap_TE1, "Overlap percentage"=(Overlap_TE1/Genes), 
                       "Analyzed in this study"=Filtered_TE1, 
                       "Overlap percentage only analyzed genes"=(Overlap_TE1/Filtered_TE1))
Data_All <- data.frame("Study"=Study, "Identified"=Genes,
                       "Overlapping"=Overlap_All, "Overlap percentage"=(Overlap_All/Genes), 
                       "Analyzed in this study"=Filtered_All, 
                       "Overlap percentage only analyzed genes"=(Overlap_All/Filtered_All))

write_xlsx(list("EE"=Data_EE, "ESR"=Data_ESR, "TE1"=Data_TE1, "All_domains"=Data_All), "./SData 9 Overlap in imprinted genes when only looking at genes present in our data.xlsx")
