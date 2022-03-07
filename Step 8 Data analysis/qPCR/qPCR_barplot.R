library(tidyverse)
library(ggpubr)
library(ggplot2)
library(readxl)
library(ggpattern)

#R version 4.1.2 (2021-11-01)
#RStudio version 2021.09.0+351 
#tidyverse version 1.3.1 
#ggplot2 version 3.3.5
#ggpubr version 0.4.0
#readxl version 1.3.1

#Import the seed phenotype frequencies of all individual siliques
All <- read_excel("./qPCR_data.xlsx", sheet="Rstudio")

All2 <- All %>% group_by(Sample) %>% summarize(Delta = mean(DeltaDeltaCt), sd=sd(DeltaDeltaCt), n = n(), se = sd / sqrt(n))

All2$Sample <- factor(All2$Sample , levels=c("ESR+", "Endosperm-", "DAL+", "All-"))
All2$GFP <- c("Neg","Pos","Neg","Pos")
All2$Domain <- c("DAL","DAL","ESR","ESR")
All2$Sample <- NULL

All2$Domain <- factor(All2$Domain , levels=c("ESR", "DAL"))
All2$GFP <- factor(All2$GFP , levels=c("Pos", "Neg"))

Plot1 <- ggplot(All2, aes(x=Domain, y=Delta, fill=GFP)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_classic() +
  scale_fill_grey() +  
  ylab("Relative GFP expression") +
  theme(text = element_text(size=10), axis.text.x = element_text(hjust=1))   +
  geom_errorbar(aes(ymin=Delta-se, ymax=Delta+se), width=.1, position=position_dodge(width=0.9)) +
  theme(legend.position = 'bottom') +
  scale_y_continuous(limits = c(0,1.25), breaks = seq(0, 1.25, 0.25))

Plot1

ggsave(filename="SFigure 1D qPCR.pdf", plot=Plot1, path='./Plots', scale = 1, height = 5, width = 5)
