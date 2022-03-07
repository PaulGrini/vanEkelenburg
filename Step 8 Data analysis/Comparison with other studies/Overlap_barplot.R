library(tidyverse)
library(ggpubr)
library(ggplot2)
library(readxl)
library(ggpattern)

#R version 4.1.1 (2021-08-10)
#RStudio version 2021.09.0+351 
#tidyverse version 1.3.1 
#ggplot2 version 3.3.5
#ggpubr version 0.4.0
#readxl version 1.3.1

# Produce in 'SData 9 Overlap in imprinted genes when only looking at genes present in our data.xlsx' a separate sheet for RStudio ('RStudio) and have it structured:
# Column A: Domains (EE, ESR, TE1, All); Column B: Study (Picard, Hornslien, Pignatta, Del_Toro); Column C: Imprinting (MEGs, PEGs); Column D: Overlap (link to the overlap percentages on each sheet);
# Column E: Type (All, Filtered). Name this 'SData 9 Overlap in imprinted genes when only looking at genes present in our data2.xlsx' (added under 'Input' on github)

#Import the seed phenotype frequencies of all individual siliques
All <- read_excel("./Input/SData 9 Overlap in imprinted genes when only looking at genes present in our data2.xlsx", sheet="RStudio")
All$Study <- factor(All$Study , levels=c("Pignatta", "Hornslien", "Del_Toro", "Picard"))
All$Domain <- factor(All$Domain , levels=c("EE", "ESR", "TE1", "All"))
MEGs <- All[which(All$Imprinting=="MEGs"),]
PEGs <- All[which(All$Imprinting=="PEGs"),]

#MEGs
Plot1 <- ggplot(MEGs, aes(x=Study, y=Overlap, fill=Type)) + 
  geom_bar(position = position_dodge(), stat = "identity") +
  facet_wrap(~Domain) +
  theme_classic() +
  scale_fill_grey(start=0.8, end=0.2) +  
  ylab("Overlapping MEGs") +
  xlab("Study") +
  theme(text = element_text(size=8), axis.text.x = element_text(hjust=1)) + 
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L), limits=c(0,0.50)) 

Plot1
ggsave(filename= "SFigure 6A Overlapping MEGs before and after filter.pdf", plot=Plot1, path='./', scale = 1, height = 7.64, width = 9.43)


#PEGs
Plot2 <- ggplot(PEGs, aes(x=Study, y=Overlap, fill=Type)) + 
  geom_bar(position = position_dodge(), stat = "identity") +
  facet_wrap(~Domain) +
  theme_classic() +
  scale_fill_grey(start=0.8, end=0.2) +  
  ylab("Overlapping PEGs") +
  xlab("Study") +
  theme(text = element_text(size=8), axis.text.x = element_text(hjust=1)) + 
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L), limits=c(0,0.60)) 

Plot2
ggsave(filename= "SFigure 6B Overlapping PEGs before and after filter.pdf", plot=Plot2, path='./', scale = 1, height = 7.64, width = 9.43)
