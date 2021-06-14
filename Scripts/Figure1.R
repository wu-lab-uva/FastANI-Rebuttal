# Generating Figure 1. The parallel function works on Mac but not on PC.
library(ape)
library(truncdist)
library(foreach)
library(doParallel)
library(data.table)
library(ggplot2)
library(ggpubr)

rm(list = ls()) #clean up the environment
registerDoParallel(cores=detectCores()) # This works for mac only
tipN = 3000 # number of tips in the unbiased tree
lambda = 19.2 # parameter for the exponential distribution
replicate = 10 # number of replicates in the simulation

# function to convert branch length into ANI
br2ani <- function (x) {    
  73.85+0.08*(100-73.85)/(x^0.63+0.08)
}

# simulation function with a tree without sampling bias
nobias_simulation <- function () {
  ani = vector()
  foreach (j=1:replicate,.combine=c) %dopar%  {
    tree = rtree(n=tipN,br=rexp(tipN,lambda))
    tipdistance = cophenetic(tree)
    tipdistance = as.vector(tipdistance[lower.tri(tipdistance,diag=FALSE)])
    ani = sapply(tipdistance,br2ani)
    return(ani)
  }
}

# simulation function with a tree sampling bias
biased_simulation <- function (N) {   # N is the number of tips that are sampled with bias
  ani = vector()
  foreach (i=1:replicate,.combine=c) %dopar% {
    tree = rtree(n=tipN,br=rexp(n=tipN,rate=lambda))
    oversamplecount = 0  #how many species will be over sampled
    while(oversamplecount < N){
      bush = rtree(2,br=rtrunc(2, spec="exp",a=0,b=qexp(0.01,rate=lambda)))
      tree = bind.tree(tree, bush, 1)
      oversamplecount= oversamplecount+1
    }
    tipdistance = cophenetic(tree)
    tipdistance = as.vector(tipdistance[lower.tri(tipdistance,diag=FALSE)])
    ani = sapply(tipdistance,br2ani)
    return(ani)
  }
}


ani_nobias <- nobias_simulation() # no tips are sampled with bias
ani_30_biased_tips <- biased_simulation(30) # 30 tips are over sampled 

hist.nobias.count = hist(ani_nobias,breaks =70:100,plot = FALSE)
hist.30_bias.count = hist(ani_30_biased_tips,breaks =70:100,plot = FALSE)
hist.nobias.data = data.frame(bias = "0",ymax= hist.nobias.count$density,ymin=0,
                              xmin=hist.nobias.count$breaks[-length(hist.nobias.count$breaks)],
                              xmax=hist.nobias.count$breaks[-1],
                              xmid = hist.nobias.count$mids)
hist.30_bias.data = data.frame(bias = "30",ymax= hist.30_bias.count$density,ymin=0,
                               xmin=hist.30_bias.count$breaks[-length(hist.30_bias.count$breaks)],
                               xmax=hist.30_bias.count$breaks[-1],
                               xmid = hist.30_bias.count$mids)
hist.all.data = rbind(hist.nobias.data,hist.30_bias.data)

bias.color = c("0"="gray","30"="gray")

hist.nobias.plot = ggplot()+
  geom_rect(mapping = aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            data = hist.nobias.data,fill="gray",col="black",alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1e-4,2.5e-5),expand = c(0,0),)+
  coord_cartesian(ylim = c(0,1e-4))+
  scale_x_continuous(limits = c(74,100),breaks=seq(70,100,5),expand = c(0,0.5))+
  scale_color_manual(values = bias.color,aesthetics = c("colour","fill"))+
  xlab(label = "ANI")+
  ylab(label = "Density")+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x =element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15))

hist.30_bias.plot = ggplot()+
  geom_rect(mapping = aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            data = hist.30_bias.data,fill="gray",col="black",alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1e-4,2.5e-5),expand = c(0,0), )+
  coord_cartesian(ylim = c(0,1e-4))+
  scale_x_continuous(limits = c(74,100),breaks=seq(70,100,5),expand = c(0,0.5))+
  scale_color_manual(values = bias.color,aesthetics = c("colour","fill"))+
  xlab(label = "ANI")+
  ylab(label = "Density")+
  #  ggtitle(label = "No bias")+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x =element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15))

plot <- function (file) { 
  dt=fread(file,header=F)
  names(dt) = c("genome1","genome2","ani")
  dt = dt[genome1 != genome2]
  plot = ggplot(dt, aes(x=ani)) + 
    geom_histogram(alpha=0.5, aes(y=..density..),color="black",fill="gray",breaks=seq(75,100,1))+
    scale_x_continuous(breaks=seq(75,100,5),expand = c(0,0.5))+
    scale_y_continuous(breaks=seq(0,2e-1,5e-2),expand = c(0,0))+
    coord_cartesian(ylim=c(0,2e-1))+
    xlab(label = "ANI")+
    ylab(label = "Density")+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text.x =element_text(face="bold", size=12),
          axis.text.y = element_text(face="bold", size=12),
          axis.title.x = element_text(face="bold", size=15),
          axis.title.y = element_text(face="bold", size=15)) 
  return(plot)
}

setwd("../Data")
random.plot = plot("ani.10genome.random2.only.10x")
representative.plot = plot("ani.10genome.representative.only")

figure1.plot = ggarrange(plotlist = list(hist.nobias.plot,hist.30_bias.plot,random.plot,representative.plot),
                         ncol = 2,nrow = 2,labels = c("a","b","c","d"),align = "hv")

figure1.plot
