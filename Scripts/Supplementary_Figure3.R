# Supplementary Figure3 - Relative tree length decay with dropping genomes per species across time

# Libraries
library(data.table)
library(tidyverse)

# For >= 10 genome
setwd('../Data')

# RTL data
out <- data.table(read.csv("SuppFig3-RTL.decay.10genomes.year.csv"))

# Change names of year
out <- data.table(out %>% mutate(year = 
                      case_when(year==2008 ~ "(1999 - 2008)",
                                year==2012 ~ "(1999 - 2012)",
                                year==2014 ~ "(1999 - 2014)",
                                year==2020 ~ "(1999 - 2018)")))

# Proportion of genomes dropped by relative branch lengths across years
ggplot() +
  geom_line(data=out, aes(x=drop, y=RTL, 
                          group=interaction(spp, year), color=as.factor(year)),
            size=1, stat="identity", alpha=0.5) +
  labs(x="Proportion of genomes dropped", y="Relative tree length") +
  facet_wrap(~year) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        strip.text = element_text(face="bold", size=12))

