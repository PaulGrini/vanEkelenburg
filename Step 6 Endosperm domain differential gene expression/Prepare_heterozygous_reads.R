library(tidyverse)
#R version 4.1.1 (2021-08-10)
#RStudio version 2021.09.0+351 
#tidyverse version 1.3.1

#Combine all EE replicates in one table
EE1 <- read_tsv("./input/count_by_gene_EE1pos.tsv",col_names = TRUE)
#Filter all reads mapped to a Col0 or Tsu1 gene
Col <- EE1 %>% filter(grepl("Col0", allele))
Tsu <- EE1 %>% filter(grepl("Tsu1", allele))
#Remove the column which specifies Col or Tsu
Col$allele <- NULL
Tsu$allele <- NULL
#Rename the column with number of read pairs to the specific ecotype
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
#Combine the two ecotypes
EE1 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
#Set al NA to 0
EE1[is.na(EE1)] <- 0
#Sum the read pairs for Col and Tsu
EE1 <- cbind(EE1, EE1pos=rowSums(EE1[,2:3]))
#Remove the separe columns for Col and Tsu
EE1$Col <- NULL
EE1$Tsu <- NULL
#Sort based on gene
EE1 <- EE1[order(EE1$gene),]

#EE3
EE3 <- read_tsv("./input/count_by_gene_EE3pos.tsv",col_names = TRUE)
Col <- EE3 %>% filter(grepl("Col0", allele))
Tsu <- EE3 %>% filter(grepl("Tsu1", allele))
Col$allele <- NULL
Tsu$allele <- NULL
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
EE3 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
EE3[is.na(EE3)] <- 0
EE3 <- cbind(EE3, EE3pos=rowSums(EE3[,2:3]))
EE3$Col <- NULL
EE3$Tsu <- NULL
EE3 <- EE3[order(EE3$gene),]

#EE4
EE4 <- read_tsv("./input/count_by_gene_EE4pos.tsv",col_names = TRUE)
Col <- EE4 %>% filter(grepl("Col0", allele))
Tsu <- EE4 %>% filter(grepl("Tsu1", allele))
Col$allele <- NULL
Tsu$allele <- NULL
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
EE4 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
EE4[is.na(EE4)] <- 0
EE4 <- cbind(EE4, EE4pos=rowSums(EE4[,2:3]))
EE4$Col <- NULL
EE4$Tsu <- NULL
EE4 <- EE4[order(EE4$gene),]

#Combine all EE samples
EE <- full_join(EE1, EE3, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., EE4, by="gene", copy=FALSE,suffix=c("1","2"))

#Combine all ESR replicates in one table
#ESR2
ESR2 <- read_tsv("./input/count_by_gene_ESR2pos.tsv",col_names = TRUE)
Col <- ESR2 %>% filter(grepl("Col0", allele))
Tsu <- ESR2 %>% filter(grepl("Tsu1", allele))
Col$allele <- NULL
Tsu$allele <- NULL
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
ESR2 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
ESR2[is.na(ESR2)] <- 0
ESR2 <- cbind(ESR2, ESR2pos=rowSums(ESR2[,2:3]))
ESR2$Col <- NULL
ESR2$Tsu <- NULL
ESR2 <- ESR2[order(ESR2$gene),]

#ESR3
ESR3 <- read_tsv("./input/count_by_gene_ESR3pos.tsv",col_names = TRUE)
Col <- ESR3 %>% filter(grepl("Col0", allele))
Tsu <- ESR3 %>% filter(grepl("Tsu1", allele))
Col$allele <- NULL
Tsu$allele <- NULL
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
ESR3 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
ESR3[is.na(ESR3)] <- 0
ESR3 <- cbind(ESR3, ESR3pos=rowSums(ESR3[,2:3]))
ESR3$Col <- NULL
ESR3$Tsu <- NULL
ESR3 <- ESR3[order(ESR3$gene),]

#ESR4
ESR4 <- read_tsv("./input/count_by_gene_ESR4pos.tsv",col_names = TRUE)
Col <- ESR4 %>% filter(grepl("Col0", allele))
Tsu <- ESR4 %>% filter(grepl("Tsu1", allele))
Col$allele <- NULL
Tsu$allele <- NULL
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
ESR4 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
ESR4[is.na(ESR4)] <- 0
ESR4 <- cbind(ESR4, ESR4pos=rowSums(ESR4[,2:3]))
ESR4$Col <- NULL
ESR4$Tsu <- NULL
ESR4 <- ESR4[order(ESR4$gene),]

#ALL ESR
ESR <- full_join(ESR2, ESR3, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., ESR4, by="gene", copy=FALSE,suffix=c("1","2"))

#Combine all DAL replicates in one table
#DAL1_5_2
DAL152 <- read_tsv("./input/count_by_gene_DAL1-5-2pos.tsv",col_names = TRUE)
Col <- DAL152 %>% filter(grepl("Col0", allele))
Tsu <- DAL152 %>% filter(grepl("Tsu1", allele))
Col$allele <- NULL
Tsu$allele <- NULL
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
DAL152 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
DAL152[is.na(DAL152)] <- 0
DAL152 <- cbind(DAL152, DAL1_5_2pos=rowSums(DAL152[,2:3]))
DAL152$Col <- NULL
DAL152$Tsu <- NULL
DAL152 <- DAL152[order(DAL152$gene),]

#DAL1_5_4
DAL154 <- read_tsv("./input/count_by_gene_DAL1-5-4pos.tsv",col_names = TRUE)
Col <- DAL154 %>% filter(grepl("Col0", allele))
Tsu <- DAL154 %>% filter(grepl("Tsu1", allele))
Col$allele <- NULL
Tsu$allele <- NULL
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
DAL154 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
DAL154[is.na(DAL154)] <- 0
DAL154 <- cbind(DAL154, DAL1_5_4pos=rowSums(DAL154[,2:3]))
DAL154$Col <- NULL
DAL154$Tsu <- NULL
DAL154 <- DAL154[order(DAL154$gene),]

#All DAL
DAL <- full_join(DAL152, DAL154, by="gene", copy=FALSE,suffix=c("1","2"))

#Combine all TE1 replicates in one table
#TE1_2
TE12 <- read_tsv("./input/count_by_gene_TE1-2pos.tsv",col_names = TRUE)
Col <- TE12 %>% filter(grepl("Col0", allele))
Tsu <- TE12 %>% filter(grepl("Tsu1", allele))
Col$allele <- NULL
Tsu$allele <- NULL
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
TE12 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
TE12[is.na(TE12)] <- 0
TE12 <- cbind(TE12, TE1_2pos=rowSums(TE12[,2:3]))
TE12$Col <- NULL
TE12$Tsu <- NULL
TE12 <- TE12[order(TE12$gene),]

#TE1_3
TE13 <- read_tsv("./input/count_by_gene_TE1-3pos.tsv",col_names = TRUE)
Col <- TE13 %>% filter(grepl("Col0", allele))
Tsu <- TE13 %>% filter(grepl("Tsu1", allele))
Col$allele <- NULL
Tsu$allele <- NULL
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
TE13 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
TE13[is.na(TE13)] <- 0
TE13 <- cbind(TE13, TE1_3pos=rowSums(TE13[,2:3]))
TE13$Col <- NULL
TE13$Tsu <- NULL
TE13 <- TE13[order(TE13$gene),]

#TE1_4
TE14 <- read_tsv("./input/count_by_gene_TE1-4pos.tsv",col_names = TRUE)
Col <- TE14 %>% filter(grepl("Col0", allele))
Tsu <- TE14 %>% filter(grepl("Tsu1", allele))
Col$allele <- NULL
Tsu$allele <- NULL
names(Col)[names(Col) == "pairs"] <- "Col"
names(Tsu)[names(Tsu) == "pairs"] <- "Tsu"
TE14 <- full_join(Col, Tsu, by="gene", copy=FALSE,suffix=c("1","2"))
TE14[is.na(TE14)] <- 0
TE14 <- cbind(TE14, TE1_4pos=rowSums(TE14[,2:3]))
TE14$Col <- NULL
TE14$Tsu <- NULL
TE14 <- TE14[order(TE14$gene),]

#ALL TE1
TE1 <- full_join(TE12, TE13, by="gene", copy=FALSE,suffix=c("1","2"))%>% full_join(., TE14, by="gene", copy=FALSE,suffix=c("1","2"))

#Exchange the NA to 0 for all merged sample tables
EE[is.na(EE)] <- 0
ESR[is.na(ESR)] <- 0
DAL[is.na(DAL)] <- 0
TE1[is.na(TE1)] <- 0

#Combine all samples in one table
All <- full_join(EE, ESR, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., DAL, by="gene", copy=FALSE,suffix=c("1","2")) %>% full_join(., TE1, by="gene", copy=FALSE,suffix=c("1","2")) 
All[is.na(All)] <- 0
All <- All[order(All$gene),]

#Mitochondrial and chloroplast genes are to be removed
All <- All[-c(24216:24424),]

#Export sorted dataset as .tsv file
write_tsv(All, "./All.tsv")
