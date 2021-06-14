# Generate Supplementary Figure 1
# FastANI and phylogenetic distance model with optimization 

# Libraries
library(ape)
library(data.table)
library(tidyverse)

rm(list=ls())
setwd("../Data")

# ANI output
ani <- readRDS("FigS1-ani.bl.new_3000_output.RDS")

# Genome trees
tree.red <- read.tree('RefSeq.3000.representative.tre')
tree.big <- read.tree('RefSeq.full.10K.tre')

# Defining variables 
true.pos <- length(ani[ani >= 95 & intra=="TRUE"]$ani) 
true.neg <- length(ani[ani < 95 & intra=="FALSE"]$ani)

false.pos <- length(ani[ani >= 95 & intra=="FALSE"]$ani) 
false.neg <- length(ani[ani < 95 & intra=="TRUE"]$ani) 

# Precision
true.pos/(true.pos + false.pos)*100

# Recall
true.pos/(true.pos + false.neg)*100

# 3K pairwise restirction
ani <- ani[g1 %in% tree.red$tip.label & g2 %in% tree.red$tip.label]

# Breaks for branch lengths
breaks <- seq(0, max(ani$value), 0.005)

# Bin distance and take average of ANI
bin.ani <- data.table(aggregate(ani$ani, 
                      by=list(cut(ani$value, 
                      breaks=seq(0,max(ani$value),0.005), 
                      labels=seq(0,max(ani$value),0.005)[1:length(breaks)-1])), median))

# rename columns
setnames(setDT(bin.ani), c("bin.dist", "avg.ani"))

# Actual data
x = seq(0, max(ani$value),0.005)[1:500]
actual.ani <- bin.ani[bin.dist%in%x]$avg.ani

# Power-law function 
predict.ani.from.branch.length = function(br,alpha,k, s){
  pred.ani = (k + ((alpha*(100-k))/((br^s) + alpha)))
  return(pred.ani)
}

# Sum of squares error generator 
exp.f <- function(params,br,ani){
  alpha=params[1]
  k=params[2]
  s=params[3]
  pred.ani = predict.ani.from.branch.length(br=br, alpha=alpha, k=k, s=s)
  sse.tot = sum((ani-pred.ani)^2)
  return(sse.tot)
}


# Optim functions
results <- optim(par=c(alpha=0.5, k=60, s=0.5), fn = exp.f, br=x, ani=actual.ani)
results$par

y = (results$par[2]+((results$par[1]*(100-results$par[2]))/((x^results$par[3]) + results$par[1])))

# ANI and branch length with model fit
ggplot() +
  geom_point(data=ani[value<2.5][!g1==g2],
             aes(x=value, y=ani), alpha=1/100) +
  geom_point(data=bin.ani[1:500], 
             aes(x=as.numeric(as.character(bin.dist)), 
                 y=avg.ani), size=1, color="#00B0F6") +
  geom_line(aes(x=x, y=y), color="#F8766D", size=1) +
  labs(x="Branch length", y="ANI") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        legend.position = "none")

