library(tidyverse)
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(writexl)
library(readxl)
library(ggpubr)
library(ggforce)

#R version 4.1.1 (2021-08-10)
#RStudio version 2021.09.0+351 
#tidyverse version 1.3.1 
#ggplot2 version 3.3.5
#VennDiagram version 1.6.20
#Writexl version 1.4.0
#dplyr version 1.0.7

#Read all tables
EE <- read.csv("../Input/pos.EE1.filtered.final.csv", header=FALSE)
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

#Sort genes on positive or negative log2FC
EE_maternal <- subset(EE, select=-c(V8,AvEExpr,t,Pvalue,B,x)) %>% .[which(.$Informative_read_log2FC > 0),] %>% .[order(.$Gene),]
EE_paternal <- subset(EE, select=-c(V8,AvEExpr,t,Pvalue,B,x)) %>% .[which(.$Informative_read_log2FC <= 0),] %>% .[order(.$Gene),]
ESR_maternal <- subset(ESR, select=-c(V8,AvESRxpr,t,Pvalue,B,x)) %>% .[which(.$Informative_read_log2FC > 0),] %>% .[order(.$Gene),]
ESR_paternal <- subset(ESR, select=-c(V8,AvESRxpr,t,Pvalue,B,x)) %>% .[which(.$Informative_read_log2FC <= 0),] %>% .[order(.$Gene),]
TE1_maternal <- subset(TE1, select=-c(V8,AvTE1xpr,t,Pvalue,B,x)) %>% .[which(.$Informative_read_log2FC > 0),] %>% .[order(.$Gene),]
TE1_paternal <- subset(TE1, select=-c(V8,AvTE1xpr,t,Pvalue,B,x)) %>% .[which(.$Informative_read_log2FC <= 0),] %>% .[order(.$Gene),]

#Sort based on gene and only keep genes that have <=30 informative reads with all 6 replicates added up
EE_maternal_not_enough_reads <- cbind(EE_maternal, All_inf_reads=rowSums(EE_maternal[,2:7])) %>% .[which(.$All_inf_reads <= 30),] %>% .[order(.$Gene),]
EE_paternal_not_enough_reads <- cbind(EE_paternal, All_inf_reads=rowSums(EE_paternal[,2:7])) %>% .[which(.$All_inf_reads <= 30),] %>% .[order(.$Gene),]
ESR_maternal_not_enough_reads <- cbind(ESR_maternal, All_inf_reads=rowSums(ESR_maternal[,2:7])) %>% .[which(.$All_inf_reads <= 30),] %>% .[order(.$Gene),]
ESR_paternal_not_enough_reads <- cbind(ESR_paternal, All_inf_reads=rowSums(ESR_paternal[,2:7])) %>% .[which(.$All_inf_reads <= 30),] %>% .[order(.$Gene),]
TE1_maternal_not_enough_reads <- cbind(TE1_maternal, All_inf_reads=rowSums(TE1_maternal[,2:7])) %>% .[which(.$All_inf_reads <= 30),] %>% .[order(.$Gene),]
TE1_paternal_not_enough_reads <- cbind(TE1_paternal, All_inf_reads=rowSums(TE1_paternal[,2:7])) %>% .[which(.$All_inf_reads <= 30),] %>% .[order(.$Gene),]

#Sort based on gene and only keep genes that have >30 informative reads with all 6 replicates added up
EE_maternal_reads_threshold <-cbind(EE_maternal, All_inf_reads=rowSums(EE_maternal[,2:7])) %>% .[which(.$All_inf_reads > 30),] %>% .[order(.$Gene),]
EE_paternal_reads_threshold <-cbind(EE_paternal, All_inf_reads=rowSums(EE_paternal[,2:7])) %>% .[which(.$All_inf_reads > 30),] %>% .[order(.$Gene),]
ESR_maternal_reads_threshold <-cbind(ESR_maternal, All_inf_reads=rowSums(ESR_maternal[,2:7])) %>% .[which(.$All_inf_reads > 30),] %>% .[order(.$Gene),]
ESR_paternal_reads_threshold <-cbind(ESR_paternal, All_inf_reads=rowSums(ESR_paternal[,2:7])) %>% .[which(.$All_inf_reads > 30),] %>% .[order(.$Gene),]
TE1_maternal_reads_threshold <-cbind(TE1_maternal, All_inf_reads=rowSums(TE1_maternal[,2:7])) %>% .[which(.$All_inf_reads > 30),] %>% .[order(.$Gene),]
TE1_paternal_reads_threshold <-cbind(TE1_paternal, All_inf_reads=rowSums(TE1_paternal[,2:7])) %>% .[which(.$All_inf_reads > 30),] %>% .[order(.$Gene),]

#Number of genes that are not parentally biased (-1 < log2 FC < 1)
EE_maternal_not_biased <- EE_maternal_reads_threshold[EE_maternal_reads_threshold$Informative_read_log2FC < 1, ] %>% .[.$Informative_read_log2FC >-1, ]  
EE_paternal_not_biased <- EE_paternal_reads_threshold[EE_paternal_reads_threshold$Informative_read_log2FC > -1, ] %>% .[.$Informative_read_log2FC >-1, ]  
ESR_maternal_not_biased <- ESR_maternal_reads_threshold[ESR_maternal_reads_threshold$Informative_read_log2FC < 1, ] %>% .[.$Informative_read_log2FC >-1, ] 
ESR_paternal_not_biased <- ESR_paternal_reads_threshold[ESR_paternal_reads_threshold$Informative_read_log2FC > -1, ] %>% .[.$Informative_read_log2FC >-1, ]  
TE1_maternal_not_biased <- TE1_maternal_reads_threshold[TE1_maternal_reads_threshold$Informative_read_log2FC < 1, ] %>% .[.$Informative_read_log2FC >-1, ]  
TE1_paternal_not_biased <- TE1_paternal_reads_threshold[TE1_paternal_reads_threshold$Informative_read_log2FC > -1, ] %>% .[.$Informative_read_log2FC >-1, ]  

#Number of genes that are parentally biased (log2 FC >= 1 or <= -1), but not significant (adj. pvalue > 0.05) 
EE_maternal_biased_not_significant <- EE_maternal_reads_threshold[EE_maternal_reads_threshold$Informative_read_log2FC >= 1, ] %>% .[which(.$Adj.Pvalue >= 0.05),] 
EE_paternal_biased_not_significant <- EE_paternal_reads_threshold[EE_paternal_reads_threshold$Informative_read_log2FC <= -1, ] %>% .[which(.$Adj.Pvalue >= 0.05),] 
ESR_maternal_biased_not_significant <- ESR_maternal_reads_threshold[ESR_maternal_reads_threshold$Informative_read_log2FC >= 1, ] %>% .[which(.$Adj.Pvalue >= 0.05),] 
ESR_paternal_biased_not_significant <- ESR_paternal_reads_threshold[ESR_paternal_reads_threshold$Informative_read_log2FC <= -1, ] %>% .[which(.$Adj.Pvalue >= 0.05),] 
TE1_maternal_biased_not_significant <- TE1_maternal_reads_threshold[TE1_maternal_reads_threshold$Informative_read_log2FC >= 1, ] %>% .[which(.$Adj.Pvalue >= 0.05),] 
TE1_paternal_biased_not_significant <- TE1_paternal_reads_threshold[TE1_paternal_reads_threshold$Informative_read_log2FC <= -1, ] %>% .[which(.$Adj.Pvalue >= 0.05),] 

#Genes that are significant imprinted
EE_significant_MEG <- EE_maternal_reads_threshold[EE_maternal_reads_threshold$Informative_read_log2FC >= 1, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
EE_significant_PEG <- EE_paternal_reads_threshold[EE_paternal_reads_threshold$Informative_read_log2FC <= -1, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
ESR_significant_MEG <- ESR_maternal_reads_threshold[ESR_maternal_reads_threshold$Informative_read_log2FC >= 1, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
ESR_significant_PEG <- ESR_paternal_reads_threshold[ESR_paternal_reads_threshold$Informative_read_log2FC <= -1, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
TE1_significant_MEG <- TE1_maternal_reads_threshold[TE1_maternal_reads_threshold$Informative_read_log2FC >= 1, ] %>% .[which(.$Adj.Pvalue < 0.05),] 
TE1_significant_PEG <- TE1_paternal_reads_threshold[TE1_paternal_reads_threshold$Informative_read_log2FC <= -1, ] %>% .[which(.$Adj.Pvalue < 0.05),] 

rpie = 1
rlabel = 1.05 * rpie # now we place labels outside of the pies

EE_pie_chart <- data.frame(Parent=c("MEG","MEG","MEG","MEG", "PEG","PEG","PEG","PEG"),
  group=c("Not enough informative reads", "Not biased", "No significant bias", "Imprinted",
          "Not enough informative reads", "Not biased", "No significant bias", "Imprinted"),
  value=c(nrow(EE_maternal_not_enough_reads), nrow(EE_maternal_not_biased), nrow(EE_maternal_biased_not_significant), nrow(EE_significant_MEG),
          nrow(EE_paternal_not_enough_reads), nrow(EE_paternal_not_biased), nrow(EE_paternal_biased_not_significant), nrow(EE_significant_PEG)))

EE_pie_chart$group <- factor(EE_pie_chart$group , levels=c("Not enough informative reads", "Not biased", "No significant bias", "Imprinted"))
EE_pies <- left_join(EE_pie_chart,
                      EE_pie_chart %>% 
                        group_by(Parent) %>%
                        summarize(Cnt_total = sum(value))) %>%
  group_by(Parent) %>%
  mutate(end_angle = 2*pi*cumsum(value)/Cnt_total,      # ending angle for each pie slice
         start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
         mid_angle = 0.5*(start_angle + end_angle))   # middle of each pie slice, for the text label

EE_pies4 <- mutate(EE_pies,
                    hjust = ifelse(mid_angle>pi, 1, 0),
                    vjust = ifelse(mid_angle<pi/2 | mid_angle>3*pi/2, 0, 1))

EE_pie <-
  ggplot(EE_pies4) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie,
                   start = start_angle, end = end_angle, fill = group))+
  geom_text(aes(x = rlabel*sin(mid_angle), y = rlabel*cos(mid_angle), label = paste0 (value),
                hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1.3, 1.3), name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.1, 1.1), name = "", breaks = NULL, labels = NULL) +
  facet_wrap(.~Parent)+
  theme_bw() +
  theme(legend.position="bottom", legend.direction="horizontal", legend.margin = margin(30, 0, 10, 0))+
  guides (fill =  guide_legend (title.theme = element_text (face = "bold")))
EE_pie

ESR_pie_chart <- data.frame(Parent=c("MEG","MEG","MEG","MEG", "PEG","PEG","PEG","PEG"),
                           group=c("Not enough informative reads", "Not biased", "No significant bias", "Imprinted",
                                   "Not enough informative reads", "Not biased", "No significant bias", "Imprinted"),
                           value=c(nrow(ESR_maternal_not_enough_reads), nrow(ESR_maternal_not_biased), nrow(ESR_maternal_biased_not_significant), nrow(ESR_significant_MEG),
                                   nrow(ESR_paternal_not_enough_reads), nrow(ESR_paternal_not_biased), nrow(ESR_paternal_biased_not_significant), nrow(ESR_significant_PEG)))

ESR_pie_chart$group <- factor(ESR_pie_chart$group , levels=c("Not enough informative reads", "Not biased", "No significant bias", "Imprinted"))
ESR_pies <- left_join(ESR_pie_chart,
                     ESR_pie_chart %>% 
                       group_by(Parent) %>%
                       summarize(Cnt_total = sum(value))) %>%
  group_by(Parent) %>%
  mutate(end_angle = 2*pi*cumsum(value)/Cnt_total,      # ending angle for each pie slice
         start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
         mid_angle = 0.5*(start_angle + end_angle))   # middle of each pie slice, for the text label

ESR_pies4 <- mutate(ESR_pies,
                   hjust = ifelse(mid_angle>pi, 1, 0),
                   vjust = ifelse(mid_angle<pi/2 | mid_angle>3*pi/2, 0, 1))

ESR_pie <-
  ggplot(ESR_pies4) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie,
                   start = start_angle, end = end_angle, fill = group))+
  geom_text(aes(x = rlabel*sin(mid_angle), y = rlabel*cos(mid_angle), label = paste0 (value),
                hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1.3, 1.3), name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.1, 1.1), name = "", breaks = NULL, labels = NULL) +
  facet_wrap(.~Parent)+
  theme_bw() +
  theme(legend.position="bottom", legend.direction="horizontal", legend.margin = margin(30, 0, 10, 0))+
  guides (fill =  guide_legend (title.theme = element_text (face = "bold")))
ESR_pie

TE1_pie_chart <- data.frame(Parent=c("MEG","MEG","MEG","MEG", "PEG","PEG","PEG","PEG"),
                           group=c("Not enough informative reads", "Not biased", "No significant bias", "Imprinted",
                                   "Not enough informative reads", "Not biased", "No significant bias", "Imprinted"),
                           value=c(nrow(TE1_maternal_not_enough_reads), nrow(TE1_maternal_not_biased), nrow(TE1_maternal_biased_not_significant), nrow(TE1_significant_MEG),
                                   nrow(TE1_paternal_not_enough_reads), nrow(TE1_paternal_not_biased), nrow(TE1_paternal_biased_not_significant), nrow(TE1_significant_PEG)))

TE1_pie_chart$group <- factor(TE1_pie_chart$group , levels=c("Not enough informative reads", "Not biased", "No significant bias", "Imprinted"))
TE1_pies <- left_join(TE1_pie_chart,
                     TE1_pie_chart %>% 
                       group_by(Parent) %>%
                       summarize(Cnt_total = sum(value))) %>%
  group_by(Parent) %>%
  mutate(end_angle = 2*pi*cumsum(value)/Cnt_total,      # ending angle for each pie slice
         start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
         mid_angle = 0.5*(start_angle + end_angle))   # middle of each pie slice, for the text label

TE1_pies4 <- mutate(TE1_pies,
                   hjust = ifelse(mid_angle>pi, 1, 0),
                   vjust = ifelse(mid_angle<pi/2 | mid_angle>3*pi/2, 0, 1))

TE1_pie <-
  ggplot(TE1_pies4, alpha=group) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie,
                   start = start_angle, end = end_angle, fill = group))+
  geom_text(aes(x = rlabel*sin(mid_angle), y = rlabel*cos(mid_angle), label = paste0 (value),
                hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1.3, 1.3), name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.1, 1.1), name = "", breaks = NULL, labels = NULL) +
  facet_wrap(.~Parent)+
  theme_bw() +
  theme(legend.position="bottom", legend.direction="horizontal", legend.margin = margin(30, 0, 10, 0))+
  guides (fill =  guide_legend (title.theme = element_text (face = "bold"))) 
TE1_pie

All_pie <- ggarrange(
  EE_pie, ESR_pie, TE1_pie, ncol=1,
  common.legend = TRUE, legend = "bottom"
)

ggsave('Pie chart EE.pdf', width=10.82, height=6.64, plot=EE_pie, device = "pdf", path="./", dpi=600)
ggsave('Pie chart ESR.pdf', width=10.82, height=6.64, plot=ESR_pie, device = "pdf", path="./", dpi=600)
ggsave('Pie chart TE1.pdf', width=10.82, height=6.64, plot=TE1_pie, device = "pdf", path="./", dpi=600)
ggsave('SFigure 5 Pie charts.pdf', width=7.82, height=11.34, plot=All_pie, device = "pdf", path="./", dpi=600)

