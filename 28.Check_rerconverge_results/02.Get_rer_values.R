library(RERconverge)

setwd("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/trees_split")
Trees=readTrees(paste('combined.tree',sep=""))
mamRERw=RERconverge::getAllResiduals(Trees, transform = "sqrt",  n.pcs = 0, use.weights = T,
                                     weights=NULL,norm="scale")

write.csv(mamRERw, 'C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/28.Check_rerconverge_results/rer_combined.csv')

setwd("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/28.Check_rerconverge_results")
pdf("rer_tree_period.pdf", height = 10, width = 10)
treePlotRers(treesObj=Trees, rermat=mamRERw, index="cwo.fas.best.fas.ClipKIT",
             type="c", nlevels=9, figwid=1.2)
dev.off()


pdf("test.pdf", height=10, width=10)
try({
  treePlotRers(treesObj=Trees, rermat=mamRERw, index="period.fas.best.fas.ClipKIT",
               type="c", nlevels=9, figwid=1.2)
})
dev.off()