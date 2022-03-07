library(tidyverse)
library(readxl)
library(gtools)
library(writexl)
#R version 4.1.1 (2021-08-10)
#RStudio version 1.4.1106 
#tidyverse version 1.3.1
#readxl version 1.3.1
#gtools version 3.9.2

#Read all tables
EE <- read.csv("./Input/pos.EE1.filtered.final.csv", header=FALSE)
names(EE)[names(EE) == "V1"] <- "Gene"
names(EE)[names(EE) == "V2"] <- "Col1"
names(EE)[names(EE) == "V3"] <- "Col2"
names(EE)[names(EE) == "V4"] <- "Col3"
names(EE)[names(EE) == "V5"] <- "Tsu1"
names(EE)[names(EE) == "V6"] <- "Tsu2"
names(EE)[names(EE) == "V7"] <- "Tsu3"
names(EE)[names(EE) == "V9"] <- "Informative_read_log2FC"
names(EE)[names(EE) == "V10"] <- "AvEExpr"
names(EE)[names(EE) == "V11"] <- "t"
names(EE)[names(EE) == "V12"] <- "Pvalue"
names(EE)[names(EE) == "V13"] <- "Adj.Pvalue"
names(EE)[names(EE) == "V14"] <- "B"
names(EE)[names(EE) == "V15"] <- "x"

ESR <- read.csv("./Input/pos.ESR.filtered.final.csv", header=FALSE)
names(ESR)[names(ESR) == "V1"] <- "Gene"
names(ESR)[names(ESR) == "V2"] <- "Col1"
names(ESR)[names(ESR) == "V3"] <- "Col2"
names(ESR)[names(ESR) == "V4"] <- "Col3"
names(ESR)[names(ESR) == "V5"] <- "Tsu1"
names(ESR)[names(ESR) == "V6"] <- "Tsu2"
names(ESR)[names(ESR) == "V7"] <- "Tsu3"
names(ESR)[names(ESR) == "V9"] <- "Informative_read_log2FC"
names(ESR)[names(ESR) == "V10"] <- "AvESRxpr"
names(ESR)[names(ESR) == "V11"] <- "t"
names(ESR)[names(ESR) == "V12"] <- "Pvalue"
names(ESR)[names(ESR) == "V13"] <- "Adj.Pvalue"
names(ESR)[names(ESR) == "V14"] <- "B"
names(ESR)[names(ESR) == "V15"] <- "x"

TE1 <- read.csv("./Input/pos.TE1.filtered.final.csv", header=FALSE)
names(TE1)[names(TE1) == "V1"] <- "Gene"
names(TE1)[names(TE1) == "V2"] <- "Col1"
names(TE1)[names(TE1) == "V3"] <- "Col2"
names(TE1)[names(TE1) == "V4"] <- "Col3"
names(TE1)[names(TE1) == "V5"] <- "Tsu1"
names(TE1)[names(TE1) == "V6"] <- "Tsu2"
names(TE1)[names(TE1) == "V7"] <- "Tsu3"
names(TE1)[names(TE1) == "V9"] <- "Informative_read_log2FC"
names(TE1)[names(TE1) == "V10"] <- "AvTE1xpr"
names(TE1)[names(TE1) == "V11"] <- "t"
names(TE1)[names(TE1) == "V12"] <- "Pvalue"
names(TE1)[names(TE1) == "V13"] <- "Adj.Pvalue"
names(TE1)[names(TE1) == "V14"] <- "B"
names(TE1)[names(TE1) == "V15"] <- "x"

#Sort based on gene and only keep genes that have >30 informative reads with all 6 replicates added up
EE <- subset(EE, select=-c(V8,AvEExpr,t,Pvalue,B,x)) %>% cbind(., All_inf_reads=rowSums(.[,2:7])) %>% .[order(.$Gene),]
ESR <- subset(ESR, select=-c(V8,AvESRxpr,t,Pvalue,B,x)) %>% cbind(., All_inf_reads=rowSums(.[,2:7])) %>% .[order(.$Gene),]
TE1 <- subset(TE1, select=-c(V8,AvTE1xpr,t,Pvalue,B,x)) %>% cbind(., All_inf_reads=rowSums(.[,2:7])) %>% .[order(.$Gene),]

#Add gene description
Description <- read_excel("./Input/All_gene_descriptions.xlsx")
Description_EE <- subset(Description, Gene %in% EE$Gene)
Description_ESR <- subset(Description, Gene %in% ESR$Gene)
Description_TE1 <- subset(Description, Gene %in% TE1$Gene)

EE <- full_join(Description_EE, EE, by="Gene", copy=FALSE,suffix=c("1","2")) 
ESR <- full_join(Description_ESR, ESR, by="Gene", copy=FALSE,suffix=c("1","2"))
TE1 <- full_join(Description_TE1, TE1, by="Gene", copy=FALSE,suffix=c("1","2"))

#Write limma output for each domain 
write_xlsx(list("EE"=EE,"ESR"=ESR, "TE1"=)), "./Output/SData 5 Limma output EE, ESR and TE1.xlsx") 

EE <- EE[which(EE$All_inf_reads > 30),]
ESR <- ESR[which(ESR$All_inf_reads > 30),]
TE1 <- TE1[which(TE1$All_inf_reads > 30),]

# Filter data for at least 31 informative read across all replicates 
write_csv(EE, "./Total inf read and ecotype specific gene filter/EE_inf_reads_and_gene_description.csv")
write_csv(ESR, "./Total inf read and ecotype specific gene filter/ESR_inf_reads_and_gene_description.csv")
write_csv(TE1, "./Total inf read and ecotype specific gene filter/TE1_inf_reads_and_gene_description.csv")
