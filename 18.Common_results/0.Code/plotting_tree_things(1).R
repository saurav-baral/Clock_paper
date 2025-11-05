
# 
setwd("G:/My Drive/Circadian Rhythm Genes Project/18.Common_results/2.State Tree")

# for this I am using phytools
# http://blog.phytools.org/2024/03/function-for-plotting-discrete-andor.html

library(phytools)

library(ape); library(phytools); library(TreeTools)
packageVersion("phytools")


tree<-read.tree("Newick Export_megpap_root_fig.nwk")
# x<-as.matrix(read.csv("x.csv",row.names=1))[,1]

if (is.null(tree$edge.length)) tree$edge.length <- rep(0.001, nrow(tree$edge))

tree <- compute.brlen(tree, method="Grafen")

pdf("State_tree_color.pdf", height = 20, width = 10)
state.data<-read.csv("two_state_species.csv",row.names=1)
fmode<-as.factor(setNames(state.data[,1],rownames(state.data)))
dotTree(tree,fmode,colors=setNames(c("blue","red"),
                                       c("diapause","no_diapause")),ftype="i",fsize=1)
dev.off()





