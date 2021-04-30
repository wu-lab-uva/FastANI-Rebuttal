# Figure2 - Relative tree length decay with dropping genomes per species

# Libraries
library(data.table)
library(tidyverse)
library(cowplot)

# For >= 10 genomes
setwd("../Data")

# RTL data
out <- data.table(read.csv("Figure2-RTL.decay.10genomes.csv"))
out[, rtl.round:= round(RTL, 2)]

# Genomes dropped at RTL = 0.95
quantile(c(out[genome.count >= 100 & rtl.round == 0.95, 
             mean(drop), spp]$V1, 0.830303030, 0.9771429))

mean(c(out[genome.count >= 100 & rtl.round == 0.95, 
               mean(drop), spp]$V1, 0.830303030, 0.9771429))

# 5 largest species
out.5 <- out[spp%in%c("Escherichia_coli", "Klebsiella_pneumoniae", 
                      "Bordetella_pertussis", "Staphylococcus_aureus", 
                      "Salmonella_enterica")]

out.5 <- data.table(out.5 %>% group_by(spp) %>% arrange(desc(genome.count)))
out <- data.table(out %>% group_by(spp) %>% arrange(desc(genome.count)))

# Proportion of genomes dropped by relative branch lengths
a <- ggplot() + 
  geom_line(data=out, aes(x=drop, y=RTL, group=spp),
            size=1, alpha=0.2, stat="identity") +
  geom_line(data=out.5, aes(x=drop, y=RTL, group=spp, color=spp), 
            size=1.1, stat = "identity") +
  scale_color_discrete(name = "Top 5 species", 
    labels=c("Escherichia_coli"="Escherichia coli", 
             "Salmonella_enterica"="Salmonella enterica",
             "Bordetella_pertussis"="Bordetella pertussis", 
             "Staphylococcus_aureus"="Staphylococcus aureus",
             "Klebsiella_pneumoniae" = "Klebsiella pneumoniae"),
    breaks=c("Escherichia_coli", "Salmonella_enterica", 
             "Bordetella_pertussis", "Staphylococcus_aureus",
             "Klebsiella_pneumoniae"))  +
  labs(x="Proportion of genomes dropped", y="Relative tree length") +
  theme_classic() +
  theme(legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold", size=13),
        legend.position = c(0.6, 0.48), 
        legend.background = element_blank(),
        axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15))

# Taxonomic information for 10K dataset
spp <- data.table(read.csv("Figure2-species.information.genome.csv"))

# Counts of genomes by species
big <- data.table(spp %>% group_by(species, taxid) %>% 
                          dplyr::summarize(counts = n()))

# Color species
big <- data.table(big %>% mutate(color=case_when(
  species=="Escherichia coli" ~ 1,
  species=="Salmonella enterica" ~ 2,
  species=="Bordetella pertussis" ~ 3,
  species=="Staphylococcus aureus" ~ 4,
  species=="Klebsiella pneumoniae" ~ 5,
  !species %in% c("Escherichia coli", "Salmonella enterica", 
                  "Bordetella pertussis", "Staphylococcus aureus", 
                  "Klebsiella pneumoniae") ~ 6)))

# Remove 2 species with no RTL decay
big <- big[!species %in% c("Rickettsia japonica", "Chlamydia muridarum")]

# Control color for columns  
col <- c("1"="#A3A500", "2"="#00B0F6", "3"="#F8766D", 
         "4"="#E76BF3", "5"="#00BF7D", "6"="grey60")

# Rank order based on genome counts 
big$taxid <- factor(big$taxid, levels=data.table(big %>% group_by(taxid) %>% 
                                               arrange(desc(counts)))$taxid)

# Genome counts for each species
b <- ggplot() +
  geom_col(data=big[counts >= 10], 
           aes(x=taxid, y=log10(counts), fill=as.factor(color)), color="black") + 
  labs(x= "Species", y="Genome counts") +
  scale_fill_manual(values = col) +
  scale_y_continuous(breaks=c(0,1,log10(25),log10(50),2,
      log10(200),log10(300),log10(400),log10(500),log10(600)),
      labels=c(0,10,25,50,100,200,300,400,500,600)) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=12, vjust=6),
        axis.title.y = element_text(face="bold", size=12),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Combine inset graphs
ggdraw() +
  draw_plot(a) +
  draw_plot(b, x = 0.07, y = 0.07, width = 0.7, height = 0.7)

