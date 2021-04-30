## Generating Supplementary Figure 2 from raw data.
library(ape)

setwd("../Data")

# read in blastn results
blast.result = read.table("E.coli.good.match.blast.txt",stringsAsFactors = FALSE,header = TRUE)
RefSeq.matched.acc = unique(blast.result$Subject_E_coli_index)
# read in GreenGene OTU map
GG.otu99.map = readRDS("GG.OTU99.map.list.RDS")
#
E.coli.otu99.acc = read.table("GreenGene.13.8.OTU99.E.coli.txt")[,1]
E.coli.otu99.list = GG.otu99.map[unique(unlist(lapply(E.coli.otu99.acc,function(x){
  which(sapply(GG.otu99.map,function(y){
  any(y==x)
}))})))]
E.coli.otu99.rep = sapply(E.coli.otu99.list,intersect,E.coli.otu99.acc)
#
RefSeq.matched.otu = unique(unlist(lapply(RefSeq.matched.acc,function(x){
  E.coli.otu99.rep[which(sapply(E.coli.otu99.list,function(y){
    any(x==y)
  }))]
})))
#
GG.full.tree = read.tree("GG.E.Coli.tre")
GG.matched.tree = drop.tip(GG.full.tree,setdiff(GG.full.tree$tip.label,RefSeq.matched.otu))
#write.tree(GG.matched.tree,"GG.E.coli.matched.tre")
# Create a pseudo tree for coloring
GG.pseudo.tree = GG.full.tree
GG.pseudo.tree$edge.length = sqrt(GG.full.tree$edge.length)
GG.matched.status = rep(0,Ntip(GG.full.tree))
GG.matched.status[match(GG.matched.tree$tip.label,GG.full.tree$tip.label)] = 1
GG.matched.status = c(GG.matched.status,
                      ace(x = GG.matched.status,phy = GG.full.tree,type = "continuous",method = "pic")$ace)
GG.edge.color = rep("red",length(GG.full.tree$edge.length))
GG.edge.color[GG.matched.status[GG.full.tree$edge[,2]]>0] = "black"
GG.tip.color = rep("red",Ntip(GG.full.tree))
GG.tip.color[match(GG.matched.tree$tip.label,GG.full.tree$tip.label)] = "black" 
##
GG.phylo.plot = plot.phylo(x = GG.pseudo.tree,type = "fan",align.tip.label = TRUE,edge.color = GG.edge.color,tip.color = GG.tip.color)
##
cover.ratio = sum(GG.matched.tree$edge.length)/sum(GG.full.tree$edge.length)
##

