# Generate Supplementary Figure 5

library(ape)
library(data.table)
library(ggplot2)

rm(list = ls()) #clean up the environment
lambda = 19.2 # parameter for the exponential distribution
replicate = 10 # number of replicates in the simulation

# function to convert branch length into ANI
br2ani <- function (x) {    
  73.85+0.08*(100-73.85)/(x^0.63+0.08)
}

# simulation function with a tree without sampling bias
simulation <- function (x) {
  ani = vector()
  for (j in 1:replicate)  {
    tree = rtree(n=x,br=rexp(x,lambda))
    tipdistance = cophenetic(tree)
    tipdistance = as.vector(tipdistance[lower.tri(tipdistance,diag=FALSE)])
    ani = append (ani, sapply(tipdistance,br2ani))
  }
  percentage = sum(ani>83 & ani<95)/length(ani)*100
  return(percentage)
}

tipN= seq(50,3000,50);
percentage = sapply(tipN, simulation)
dt = data.table(tipN,percentage)

ggplot(dt,aes(tipN,percentage))+
  geom_point(shape = 21,fill="red",color="black")+
  xlab(label = "Number of genomes")+
  ylab(label = "% of ANI in [83%-95%]")+
  scale_x_continuous(limits = c(0,3100),breaks=seq(0,3000,500),expand = c(0,0))+
  scale_y_continuous(limits=c(0,0.55), expand=c(0,0))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x =element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15))

