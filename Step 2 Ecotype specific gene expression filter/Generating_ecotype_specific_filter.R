# Generating SData3 containing the ecotype specific filter lists and the raw data from which this was derived

library(tidyverse)
library(writexl)
#R version 4.1.1 (2021-08-10)
#RStudio version 1.4.1106 
#tidyverse version 1.3.1 
#Writexl version 1.4.0

Early <- read_tsv("./Merged_read_counts/Early_all_reads.tsv", col_names=TRUE)
names(Early)[names(Early) == "gene"] <- "Gene"
Late <- read_tsv("./Merged_read_counts/Late_all_reads.tsv", col_names=TRUE)
names(Late)[names(Late) == "gene"] <- "Gene"
Early_DESeq <- read.csv("./Output/Early_Col_vs_Tsu.csv", header=TRUE)
names(Early_DESeq)[names(Early_DESeq) == "X"] <- "Gene"
Late_DESeq <- read.csv("./Output/Late_Col_vs_Tsu.csv", header=TRUE)
names(Late_DESeq)[names(Late_DESeq) == "X"] <- "Gene"

Early_raw <- full_join(Early, Early_DESeq, by="Gene", copy=FALSE,suffix=c("1","2"))
Late_raw <- full_join(Late, Late_DESeq, by="Gene", copy=FALSE,suffix=c("1","2"))

Early_raw <- Early_raw[,c(1:7,9,12)]
Late_raw <- Late_raw[,c(1:7,9,12)]

Early <- Early_raw[which(Early_raw$padj > 0.05),]
Late <- Late_raw[which(Late_raw$padj > 0.05),]

Early <- Early[,c(1)]
Late <- Late[,c(1)]

write_tsv(Early, file = "./Output/Early_ecotype_pass_filter.tsv")
write_tsv(Late, file = "./Output/Late_ecotype_pass_filter.tsv")

write_xlsx(list("Early"=Early, "Late"=Late, "Early raw"=Early_raw, "Late raw"=Late_raw), "./Output/SData 4 Ecotype specific gene filter.xlsx")

