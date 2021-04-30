# Generate Figure 3
# Intraspecific ANI distribution for species with 100+ genomes in Jain et al, 2018.

# Libraries
library(data.table)
library(tidyverse)
library(taxize)
library(plyr)

rm(list=ls())
setwd("../Data/")
df = readRDS("Figure3.RDS")

# Color phyla
df <- data.table(df %>% mutate(color=case_when(
  phylum=="Proteobacteria" ~ "#F8766D",
  phylum=="Firmicutes" ~ "#B79F00",
  phylum=="Bacteroidetes" ~ "#00BA38",
  phylum=="Chlamydiae" ~ "#00BFC4",
  phylum=="Spirochaetes" ~ "#619CFF",
  phylum=="Actinobacteria" ~ "#F564E3")))

# Color for each species
col <- df %>% group_by(species) %>% distinct(color)

breaks <- seq(75, 100, by=5)
ggplot(df, aes(x=species,y=ani,fill=freq))+
  geom_tile()+
  coord_flip()+
  scale_fill_gradientn(colors=c("white","gray","red"), values=c(0,1e-8,1))+
  geom_hline(yintercept = 95, linetype=2)+
  geom_text(aes(y=70, label=paste(N)), size=3)+
  scale_y_continuous(breaks=c(70, breaks),)+
  labs(x="Species", y="ANI") +
  theme_classic()+
  theme(axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold.italic", size=9, color=col$color),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15)) 

ggsave("Figure3.jpeg",width=8.5, height=11, unit="in",device="jpeg")


