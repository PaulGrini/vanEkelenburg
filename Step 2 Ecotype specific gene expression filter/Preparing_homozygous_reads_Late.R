#Collect all tsv files with read counts in one folder '/Counts'
#Collect output tsv files with read counts from all replicates in '/Merged_read_counts'
#This R script merges all individual tsv files for homozygous reads at 4 DAP with read counts together

library(tidyverse)
#R version 4.1.1 (2021-08-10)
#RStudio version 1.4.1106 
#tidyverse version 1.3.1

#Import reads
#Combine all EE replicates in one table
Col0_Late_BR1_20190705 <- read_tsv("./Counts/07-Col-0-Late-BR1.20190705.tsv",col_names = TRUE)
Col0_Late_BR1_20190914 <- read_tsv("./Counts/07-Col-0-Late-BR1.20190914.tsv",col_names = TRUE)
Col0_Late_BR1_20190705$`indel-free` <- NULL
Col0_Late_BR1_20190914$`indel-free` <- NULL
Col0_Late_BR1_20190705$spliced <- NULL
Col0_Late_BR1_20190914$spliced <- NULL
#Split Col and Tsu
#Filter all reads mapped to a Col0 or Tsu1 gene
Col20190705 <- Col0_Late_BR1_20190705 %>% filter(grepl("Col0", allele))
Tsu20190705 <- Col0_Late_BR1_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Col0_Late_BR1_20190914 %>% filter(grepl("Col0", allele))
Tsu20190914 <- Col0_Late_BR1_20190914 %>% filter(grepl("Tsu1", allele))
#Remove the column which specifies Col or Tsu
Col20190705$allele <- NULL
Tsu20190705$allele <- NULL
Col20190914$allele <- NULL
Tsu20190914$allele <- NULL
#Rename the column with number of read count to the specific ecotype
names(Col20190705)[names(Col20190705) == "pairs"] <- "Col20190705"
names(Tsu20190705)[names(Tsu20190705) == "pairs"] <- "Tsu20190705"
names(Col20190914)[names(Col20190914) == "pairs"] <- "Col20190914"
names(Tsu20190914)[names(Tsu20190914) == "pairs"] <- "Tsu20190914"
#Combine the two ecotypes
ColBR1 <- full_join(Col20190705, Col20190914, by="gene", copy=FALSE,suffix=c("1","2"))
Tsu1BR1 <- full_join(Tsu20190705, Tsu20190914, by="gene", copy=FALSE,suffix=c("1","2"))
#Set al NA to 0
ColBR1[is.na(ColBR1)] <- 0
Tsu1BR1[is.na(Tsu1BR1)] <- 0
#Sum the read pairs for Col and Tsu
ColBR1 <- cbind(ColBR1, Col_BR1=rowSums(ColBR1[,2:3]))
Tsu1BR1 <- cbind(Tsu1BR1, Tsu1_BR1=rowSums(Tsu1BR1[,2:3]))
#Remove the separe columns for Col and Tsu
ColBR1$Col20190705 <- NULL
ColBR1$Col20190914 <- NULL
Tsu1BR1$Tsu20190705 <- NULL
Tsu1BR1$Tsu20190914 <- NULL
ColLateBR1 <- full_join(ColBR1, Tsu1BR1, by="gene", copy=FALSE,suffix=c("1","2"))
ColLateBR1[is.na(ColLateBR1)] <- 0
ColLateBR1 <- cbind(ColLateBR1, ColLateBR1=rowSums(ColLateBR1[,2:3]))
ColLateBR1$Col_BR1 <- NULL
ColLateBR1$Tsu1_BR1 <- NULL

#ColLateBR2
Col0_Late_BR2_20190705 <- read_tsv("./Counts/08-Col-0-Late-BR2.20190705.tsv",col_names = TRUE)
Col0_Late_BR2_20190914 <- read_tsv("./Counts/08-Col-0-Late-BR2.20190914.tsv",col_names = TRUE)
Col0_Late_BR2_20190705$`indel-free` <- NULL
Col0_Late_BR2_20190914$`indel-free` <- NULL
Col0_Late_BR2_20190705$spliced <- NULL
Col0_Late_BR2_20190914$spliced <- NULL
Col20190705 <- Col0_Late_BR2_20190705 %>% filter(grepl("Col0", allele))
Tsu20190705 <- Col0_Late_BR2_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Col0_Late_BR2_20190914 %>% filter(grepl("Col0", allele))
Tsu20190914 <- Col0_Late_BR2_20190914 %>% filter(grepl("Tsu1", allele))
Col20190705$allele <- NULL
Tsu20190705$allele <- NULL
Col20190914$allele <- NULL
Tsu20190914$allele <- NULL
names(Col20190705)[names(Col20190705) == "pairs"] <- "Col20190705"
names(Tsu20190705)[names(Tsu20190705) == "pairs"] <- "Tsu20190705"
names(Col20190914)[names(Col20190914) == "pairs"] <- "Col20190914"
names(Tsu20190914)[names(Tsu20190914) == "pairs"] <- "Tsu20190914"
ColBR2 <- full_join(Col20190705, Col20190914, by="gene", copy=FALSE,suffix=c("1","2"))
Tsu1BR2 <- full_join(Tsu20190705, Tsu20190914, by="gene", copy=FALSE,suffix=c("1","2"))
ColBR2[is.na(ColBR2)] <- 0
Tsu1BR2[is.na(Tsu1BR2)] <- 0
ColBR2 <- cbind(ColBR2, Col_BR2=rowSums(ColBR2[,2:3]))
Tsu1BR2 <- cbind(Tsu1BR2, Tsu1_BR2=rowSums(Tsu1BR2[,2:3]))
ColBR2$Col20190705 <- NULL
ColBR2$Col20190914 <- NULL
Tsu1BR2$Tsu20190705 <- NULL
Tsu1BR2$Tsu20190914 <- NULL
ColLateBR2 <- full_join(ColBR2, Tsu1BR2, by="gene", copy=FALSE,suffix=c("1","2"))
ColLateBR2[is.na(ColLateBR2)] <- 0
ColLateBR2 <- cbind(ColLateBR2, ColLateBR2=rowSums(ColLateBR2[,2:3]))
ColLateBR2$Col_BR2 <- NULL
ColLateBR2$Tsu1_BR2 <- NULL

#ColLateBR3
Col0_Late_BR3_20190705 <- read_tsv("./Counts/09-Col-0-Late-BR3.20190705.tsv",col_names = TRUE)
Col0_Late_BR3_20190914 <- read_tsv("./Counts/09-Col-0-Late-BR3.20190914.tsv",col_names = TRUE)
Col0_Late_BR3_20190705$`indel-free` <- NULL
Col0_Late_BR3_20190914$`indel-free` <- NULL
Col0_Late_BR3_20190705$spliced <- NULL
Col0_Late_BR3_20190914$spliced <- NULL
Col20190705 <- Col0_Late_BR3_20190705 %>% filter(grepl("Col0", allele))
Tsu20190705 <- Col0_Late_BR3_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Col0_Late_BR3_20190914 %>% filter(grepl("Col0", allele))
Tsu20190914 <- Col0_Late_BR3_20190914 %>% filter(grepl("Tsu1", allele))
Col20190705$allele <- NULL
Tsu20190705$allele <- NULL
Col20190914$allele <- NULL
Tsu20190914$allele <- NULL
names(Col20190705)[names(Col20190705) == "pairs"] <- "Col20190705"
names(Tsu20190705)[names(Tsu20190705) == "pairs"] <- "Tsu20190705"
names(Col20190914)[names(Col20190914) == "pairs"] <- "Col20190914"
names(Tsu20190914)[names(Tsu20190914) == "pairs"] <- "Tsu20190914"
ColBR3 <- full_join(Col20190705, Col20190914, by="gene", copy=FALSE,suffix=c("1","2"))
Tsu1BR3 <- full_join(Tsu20190705, Tsu20190914, by="gene", copy=FALSE,suffix=c("1","2"))
ColBR3[is.na(ColBR3)] <- 0
Tsu1BR3[is.na(Tsu1BR3)] <- 0
ColBR3 <- cbind(ColBR3, Col_BR3=rowSums(ColBR3[,2:3]))
Tsu1BR3 <- cbind(Tsu1BR3, Tsu1_BR3=rowSums(Tsu1BR3[,2:3]))
ColBR3$Col20190705 <- NULL
ColBR3$Col20190914 <- NULL
Tsu1BR3$Tsu20190705 <- NULL
Tsu1BR3$Tsu20190914 <- NULL
ColLateBR3 <- full_join(ColBR3, Tsu1BR3, by="gene", copy=FALSE,suffix=c("1","2"))
ColLateBR3[is.na(ColLateBR3)] <- 0
ColLateBR3 <- cbind(ColLateBR3, ColLateBR3=rowSums(ColLateBR3[,2:3]))
ColLateBR3$Col_BR3 <- NULL
ColLateBR3$Tsu1_BR3 <- NULL

#TsuLateBR1
Tsu1_Late_BR1_20190705 <- read_tsv("./Counts/10-Tsu-1-Late-BR1.20190705.tsv",col_names = TRUE)
Tsu1_Late_BR1_20190914 <- read_tsv("./Counts/10-Tsu-1-Late-BR1.20190914.tsv",col_names = TRUE)
Tsu1_Late_BR1_20190705$`indel-free` <- NULL
Tsu1_Late_BR1_20190914$`indel-free` <- NULL
Tsu1_Late_BR1_20190705$spliced <- NULL
Tsu1_Late_BR1_20190914$spliced <- NULL
Col20190705 <- Tsu1_Late_BR1_20190705 %>% filter(grepl("Tsu1", allele))
Tsu20190705 <- Tsu1_Late_BR1_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Tsu1_Late_BR1_20190914 %>% filter(grepl("Tsu1", allele))
Tsu20190914 <- Tsu1_Late_BR1_20190914 %>% filter(grepl("Tsu1", allele))
Col20190705$allele <- NULL
Tsu20190705$allele <- NULL
Col20190914$allele <- NULL
Tsu20190914$allele <- NULL
names(Col20190705)[names(Col20190705) == "pairs"] <- "Col20190705"
names(Tsu20190705)[names(Tsu20190705) == "pairs"] <- "Tsu20190705"
names(Col20190914)[names(Col20190914) == "pairs"] <- "Col20190914"
names(Tsu20190914)[names(Tsu20190914) == "pairs"] <- "Tsu20190914"
ColBR1 <- full_join(Col20190705, Col20190914, by="gene", copy=FALSE,suffix=c("1","2"))
Tsu1BR1 <- full_join(Tsu20190705, Tsu20190914, by="gene", copy=FALSE,suffix=c("1","2"))
ColBR1[is.na(ColBR1)] <- 0
Tsu1BR1[is.na(Tsu1BR1)] <- 0
ColBR1 <- cbind(ColBR1, Col_BR1=rowSums(ColBR1[,2:3]))
Tsu1BR1 <- cbind(Tsu1BR1, Tsu1_BR1=rowSums(Tsu1BR1[,2:3]))
ColBR1$Col20190705 <- NULL
ColBR1$Col20190914 <- NULL
Tsu1BR1$Tsu20190705 <- NULL
Tsu1BR1$Tsu20190914 <- NULL
TsuLateBR1 <- full_join(ColBR1, Tsu1BR1, by="gene", copy=FALSE,suffix=c("1","2"))
TsuLateBR1[is.na(TsuLateBR1)] <- 0
TsuLateBR1 <- cbind(TsuLateBR1, TsuLateBR1=rowSums(TsuLateBR1[,2:3]))
TsuLateBR1$Col_BR1 <- NULL
TsuLateBR1$Tsu1_BR1 <- NULL

#TsuLateBR2
Tsu1_Late_BR2_20190705 <- read_tsv("./Counts/11-Tsu-1-Late-BR2.20190705.tsv",col_names = TRUE)
Tsu1_Late_BR2_20190914 <- read_tsv("./Counts/11-Tsu-1-Late-BR2.20190914.tsv",col_names = TRUE)
Tsu1_Late_BR2_20190705$`indel-free` <- NULL
Tsu1_Late_BR2_20190914$`indel-free` <- NULL
Tsu1_Late_BR2_20190705$spliced <- NULL
Tsu1_Late_BR2_20190914$spliced <- NULL
Col20190705 <- Tsu1_Late_BR2_20190705 %>% filter(grepl("Tsu1", allele))
Tsu20190705 <- Tsu1_Late_BR2_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Tsu1_Late_BR2_20190914 %>% filter(grepl("Tsu1", allele))
Tsu20190914 <- Tsu1_Late_BR2_20190914 %>% filter(grepl("Tsu1", allele))
Col20190705$allele <- NULL
Tsu20190705$allele <- NULL
Col20190914$allele <- NULL
Tsu20190914$allele <- NULL
names(Col20190705)[names(Col20190705) == "pairs"] <- "Col20190705"
names(Tsu20190705)[names(Tsu20190705) == "pairs"] <- "Tsu20190705"
names(Col20190914)[names(Col20190914) == "pairs"] <- "Col20190914"
names(Tsu20190914)[names(Tsu20190914) == "pairs"] <- "Tsu20190914"
ColBR2 <- full_join(Col20190705, Col20190914, by="gene", copy=FALSE,suffix=c("1","2"))
Tsu1BR2 <- full_join(Tsu20190705, Tsu20190914, by="gene", copy=FALSE,suffix=c("1","2"))
ColBR2[is.na(ColBR2)] <- 0
Tsu1BR2[is.na(Tsu1BR2)] <- 0
ColBR2 <- cbind(ColBR2, Col_BR2=rowSums(ColBR2[,2:3]))
Tsu1BR2 <- cbind(Tsu1BR2, Tsu1_BR2=rowSums(Tsu1BR2[,2:3]))
ColBR2$Col20190705 <- NULL
ColBR2$Col20190914 <- NULL
Tsu1BR2$Tsu20190705 <- NULL
Tsu1BR2$Tsu20190914 <- NULL
TsuLateBR2 <- full_join(ColBR2, Tsu1BR2, by="gene", copy=FALSE,suffix=c("1","2"))
TsuLateBR2[is.na(TsuLateBR2)] <- 0
TsuLateBR2 <- cbind(TsuLateBR2, TsuLateBR2=rowSums(TsuLateBR2[,2:3]))
TsuLateBR2$Col_BR2 <- NULL
TsuLateBR2$Tsu1_BR2 <- NULL

#TsuLateBR3
Tsu1_Late_BR3_20190705 <- read_tsv("./Counts/12-Tsu-1-Late-BR3.20190705.tsv",col_names = TRUE)
Tsu1_Late_BR3_20190914 <- read_tsv("./Counts/12-Tsu-1-Late-BR3.20190914.tsv",col_names = TRUE)
Tsu1_Late_BR3_20190705$`indel-free` <- NULL
Tsu1_Late_BR3_20190914$`indel-free` <- NULL
Tsu1_Late_BR3_20190705$spliced <- NULL
Tsu1_Late_BR3_20190914$spliced <- NULL
Col20190705 <- Tsu1_Late_BR3_20190705 %>% filter(grepl("Tsu1", allele))
Tsu20190705 <- Tsu1_Late_BR3_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Tsu1_Late_BR3_20190914 %>% filter(grepl("Tsu1", allele))
Tsu20190914 <- Tsu1_Late_BR3_20190914 %>% filter(grepl("Tsu1", allele))
Col20190705$allele <- NULL
Tsu20190705$allele <- NULL
Col20190914$allele <- NULL
Tsu20190914$allele <- NULL
names(Col20190705)[names(Col20190705) == "pairs"] <- "Col20190705"
names(Tsu20190705)[names(Tsu20190705) == "pairs"] <- "Tsu20190705"
names(Col20190914)[names(Col20190914) == "pairs"] <- "Col20190914"
names(Tsu20190914)[names(Tsu20190914) == "pairs"] <- "Tsu20190914"
ColBR3 <- full_join(Col20190705, Col20190914, by="gene", copy=FALSE,suffix=c("1","2"))
Tsu1BR3 <- full_join(Tsu20190705, Tsu20190914, by="gene", copy=FALSE,suffix=c("1","2"))
ColBR3[is.na(ColBR3)] <- 0
Tsu1BR3[is.na(Tsu1BR3)] <- 0
ColBR3 <- cbind(ColBR3, Col_BR3=rowSums(ColBR3[,2:3]))
Tsu1BR3 <- cbind(Tsu1BR3, Tsu1_BR3=rowSums(Tsu1BR3[,2:3]))
ColBR3$Col20190705 <- NULL
ColBR3$Col20190914 <- NULL
Tsu1BR3$Tsu20190705 <- NULL
Tsu1BR3$Tsu20190914 <- NULL
TsuLateBR3 <- full_join(ColBR3, Tsu1BR3, by="gene", copy=FALSE,suffix=c("1","2"))
TsuLateBR3[is.na(TsuLateBR3)] <- 0
TsuLateBR3 <- cbind(TsuLateBR3, TsuLateBR3=rowSums(TsuLateBR3[,2:3]))
TsuLateBR3$Col_BR3 <- NULL
TsuLateBR3$Tsu1_BR3 <- NULL

Late_All <- full_join(ColLateBR1, ColLateBR2, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., ColLateBR3, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., TsuLateBR1, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., TsuLateBR2, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., TsuLateBR3, by="gene", copy=FALSE,suffix=c("1","2"))
Late_All[is.na(Late_All)] <- 0
Late_All <- Late_All[order(Late_All$gene),]

write_tsv(Late_All, "./Merged_read_counts/Late_all_reads.tsv")
