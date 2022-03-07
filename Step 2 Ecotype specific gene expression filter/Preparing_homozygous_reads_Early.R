#Collect all tsv files with read counts in one folder '/Counts'
#Collect output tsv files with read counts from all replicates in '/Merged_read_counts'
#This R script merges all individual tsv files for homozygous reads at 4 DAP with read counts together

library(tidyverse)
#R version 4.1.1 (2021-08-10)
#RStudio version 1.4.1106 
#tidyverse version 1.3.1

#Import reads
#Combine all EE replicates in one table
Col0_Early_BR1_20190705 <- read_tsv("./Counts/01-Col-0-Early-BR1.20190705.tsv",col_names = TRUE)
Col0_Early_BR1_20190914 <- read_tsv("./Counts/01-Col-0-Early-BR1.20190914.tsv",col_names = TRUE)
#Remove the column which specifies indel-free and spliced
Col0_Early_BR1_20190705$`indel-free` <- NULL
Col0_Early_BR1_20190914$`indel-free` <- NULL
Col0_Early_BR1_20190705$spliced <- NULL
Col0_Early_BR1_20190914$spliced <- NULL

#Split Col and Tsu
#Filter all reads mapped to a Col0 or Tsu1 gene
Col20190705 <- Col0_Early_BR1_20190705 %>% filter(grepl("Col0", allele))
Tsu20190705 <- Col0_Early_BR1_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Col0_Early_BR1_20190914 %>% filter(grepl("Col0", allele))
Tsu20190914 <- Col0_Early_BR1_20190914 %>% filter(grepl("Tsu1", allele))
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
ColEarlyBR1 <- full_join(ColBR1, Tsu1BR1, by="gene", copy=FALSE,suffix=c("1","2"))
ColEarlyBR1[is.na(ColEarlyBR1)] <- 0
ColEarlyBR1 <- cbind(ColEarlyBR1, ColEarlyBR1=rowSums(ColEarlyBR1[,2:3]))
ColEarlyBR1$Col_BR1 <- NULL
ColEarlyBR1$Tsu1_BR1 <- NULL

#ColEarlyBR2
Col0_Early_BR2_20190705 <- read_tsv("./Counts/02-Col-0-Early-BR2.20190705.tsv",col_names = TRUE)
Col0_Early_BR2_20190914 <- read_tsv("./Counts/02-Col-0-Early-BR2.20190914.tsv",col_names = TRUE)
#Remove the column which specifies indel-free and spliced
Col0_Early_BR2_20190705$`indel-free` <- NULL
Col0_Early_BR2_20190914$`indel-free` <- NULL
Col0_Early_BR2_20190705$spliced <- NULL
Col0_Early_BR2_20190914$spliced <- NULL
Col20190705 <- Col0_Early_BR2_20190705 %>% filter(grepl("Col0", allele))
Tsu20190705 <- Col0_Early_BR2_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Col0_Early_BR2_20190914 %>% filter(grepl("Col0", allele))
Tsu20190914 <- Col0_Early_BR2_20190914 %>% filter(grepl("Tsu1", allele))
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
ColEarlyBR2 <- full_join(ColBR2, Tsu1BR2, by="gene", copy=FALSE,suffix=c("1","2"))
ColEarlyBR2[is.na(ColEarlyBR2)] <- 0
ColEarlyBR2 <- cbind(ColEarlyBR2, ColEarlyBR2=rowSums(ColEarlyBR2[,2:3]))
ColEarlyBR2$Col_BR2 <- NULL
ColEarlyBR2$Tsu1_BR2 <- NULL

#ColEarlyBR3
Col0_Early_BR3_20190705 <- read_tsv("./Counts/03-Col-0-Early-BR3.20190705.tsv",col_names = TRUE)
Col0_Early_BR3_20190914 <- read_tsv("./Counts/03-Col-0-Early-BR3.20190914.tsv",col_names = TRUE)
#Remove the column which specifies indel-free and spliced
Col0_Early_BR3_20190705$`indel-free` <- NULL
Col0_Early_BR3_20190914$`indel-free` <- NULL
Col0_Early_BR3_20190705$spliced <- NULL
Col0_Early_BR3_20190914$spliced <- NULL
Col20190705 <- Col0_Early_BR3_20190705 %>% filter(grepl("Col0", allele))
Tsu20190705 <- Col0_Early_BR3_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Col0_Early_BR3_20190914 %>% filter(grepl("Col0", allele))
Tsu20190914 <- Col0_Early_BR3_20190914 %>% filter(grepl("Tsu1", allele))
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
ColEarlyBR3 <- full_join(ColBR3, Tsu1BR3, by="gene", copy=FALSE,suffix=c("1","2"))
ColEarlyBR3[is.na(ColEarlyBR3)] <- 0
ColEarlyBR3 <- cbind(ColEarlyBR3, ColEarlyBR3=rowSums(ColEarlyBR3[,2:3]))
ColEarlyBR3$Col_BR3 <- NULL
ColEarlyBR3$Tsu1_BR3 <- NULL

#TsuEarlyBR1
Tsu1_Early_BR1_20190705 <- read_tsv("./Counts/04-Tsu-1-Early-BR1.20190705.tsv",col_names = TRUE)
Tsu1_Early_BR1_20190914 <- read_tsv("./Counts/04-Tsu-1-Early-BR1.20190914.tsv",col_names = TRUE)
#Remove the column which specifies indel-free and spliced
Tsu1_Early_BR1_20190705$`indel-free` <- NULL
Tsu1_Early_BR1_20190914$`indel-free` <- NULL
Tsu1_Early_BR1_20190705$spliced <- NULL
Tsu1_Early_BR1_20190914$spliced <- NULL
Col20190705 <- Tsu1_Early_BR1_20190705 %>% filter(grepl("Tsu1", allele))
Tsu20190705 <- Tsu1_Early_BR1_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Tsu1_Early_BR1_20190914 %>% filter(grepl("Tsu1", allele))
Tsu20190914 <- Tsu1_Early_BR1_20190914 %>% filter(grepl("Tsu1", allele))
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
TsuEarlyBR1 <- full_join(ColBR1, Tsu1BR1, by="gene", copy=FALSE,suffix=c("1","2"))
TsuEarlyBR1[is.na(TsuEarlyBR1)] <- 0
TsuEarlyBR1 <- cbind(TsuEarlyBR1, TsuEarlyBR1=rowSums(TsuEarlyBR1[,2:3]))
TsuEarlyBR1$Col_BR1 <- NULL
TsuEarlyBR1$Tsu1_BR1 <- NULL

#TsuEarlyBR2
Tsu1_Early_BR2_20190705 <- read_tsv("./Counts/05-Tsu-1-Early-BR2.20190705.tsv",col_names = TRUE)
Tsu1_Early_BR2_20190914 <- read_tsv("./Counts/05-Tsu-1-Early-BR2.20190914.tsv",col_names = TRUE)
Tsu1_Early_BR2_20190705$`indel-free` <- NULL
Tsu1_Early_BR2_20190914$`indel-free` <- NULL
Tsu1_Early_BR2_20190705$spliced <- NULL
Tsu1_Early_BR2_20190914$spliced <- NULL
Col20190705 <- Tsu1_Early_BR2_20190705 %>% filter(grepl("Tsu1", allele))
Tsu20190705 <- Tsu1_Early_BR2_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Tsu1_Early_BR2_20190914 %>% filter(grepl("Tsu1", allele))
Tsu20190914 <- Tsu1_Early_BR2_20190914 %>% filter(grepl("Tsu1", allele))
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
TsuEarlyBR2 <- full_join(ColBR2, Tsu1BR2, by="gene", copy=FALSE,suffix=c("1","2"))
TsuEarlyBR2[is.na(TsuEarlyBR2)] <- 0
TsuEarlyBR2 <- cbind(TsuEarlyBR2, TsuEarlyBR2=rowSums(TsuEarlyBR2[,2:3]))
TsuEarlyBR2$Col_BR2 <- NULL
TsuEarlyBR2$Tsu1_BR2 <- NULL

#TsuEarlyBR3
Tsu1_Early_BR3_20190705 <- read_tsv("./Counts/06-Tsu-1-Early-BR3.20190705.tsv",col_names = TRUE)
Tsu1_Early_BR3_20190914 <- read_tsv("./Counts/06-Tsu-1-Early-BR3.20190914.tsv",col_names = TRUE)
Tsu1_Early_BR3_20190705$`indel-free` <- NULL
Tsu1_Early_BR3_20190914$`indel-free` <- NULL
Tsu1_Early_BR3_20190705$spliced <- NULL
Tsu1_Early_BR3_20190914$spliced <- NULL
Col20190705 <- Tsu1_Early_BR3_20190705 %>% filter(grepl("Tsu1", allele))
Tsu20190705 <- Tsu1_Early_BR3_20190705 %>% filter(grepl("Tsu1", allele))
Col20190914 <- Tsu1_Early_BR3_20190914 %>% filter(grepl("Tsu1", allele))
Tsu20190914 <- Tsu1_Early_BR3_20190914 %>% filter(grepl("Tsu1", allele))
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
TsuEarlyBR3 <- full_join(ColBR3, Tsu1BR3, by="gene", copy=FALSE,suffix=c("1","2"))
TsuEarlyBR3[is.na(TsuEarlyBR3)] <- 0
TsuEarlyBR3 <- cbind(TsuEarlyBR3, TsuEarlyBR3=rowSums(TsuEarlyBR3[,2:3]))
TsuEarlyBR3$Col_BR3 <- NULL
TsuEarlyBR3$Tsu1_BR3 <- NULL

Early_All <- full_join(ColEarlyBR1, ColEarlyBR2, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., ColEarlyBR3, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., TsuEarlyBR1, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., TsuEarlyBR2, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., TsuEarlyBR3, by="gene", copy=FALSE,suffix=c("1","2"))
Early_All[is.na(Early_All)] <- 0
Early_All <- Early_All[order(Early_All$gene),]

write_tsv(Early_All, "./Merged_read_counts/Early_all_reads.tsv")
