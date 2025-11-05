library(RERconverge)
setwd("C:/Users/Saurav Baral/Desktop/r_related")
toyTrees=readTrees(paste('new_domain_2.tre',sep=""), max.read = 1800)

toyTrees=readTrees(paste("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/5.From_Busco/combined_tree_old.txt",sep=""), max.read = 200)
new_toy_tree  =fixPseudoroot(toyTrees)

wtiplabel <- which(toyTrees$masterTree$tip.label=="tiplabel.manuallyselectedrootspecies")
newMaster <- unroot(reroot(toyTrees$masterTree,c("Calephelis_nemesis")))

mamRERw=RERconverge::getAllResiduals(toyTrees,
                                     transform = "sqrt", weighted = T, scale = T)
if (!exists('plotTreeHighlightBranches')) {
  repopath = '~/repos/RERconverge' #modify for other machines
  suppressMessages(source(paste(repopath,'/R/plottingFuncs.R',sep='')))
}

noneutherians <- c("Nymphalis_c","Nymphalis_polychloros","Nymphalis_io","Nymphalis_urticae")
par(mfrow=c(1,2))
avgtree=plotTreeHighlightBranches(toyTrees$masterTree,
                                  hlspecies=c("Napeogenes_inachia","Pararge_aegeria"), hlcols=c("blue","red"),
                                  main="Average tree") #plot average tree
bend3tree=plotTreeHighlightBranches(toyTrees$trees$period,
                                    hlspecies=c("Napeogenes_inachia","Pararge_aegeria"), hlcols=c("blue","red"),
                                    main="BEND3 tree") #plot individual gene tree
par(mfrow=c(1,1))
phenvExample <- foreground2Paths(c("Napeogenes_inachia","Pararge_aegeria"),toyTrees,clade="terminal")
plotRers(mamRERw,"period",phenv=phenvExample)



marineb=read.tree(paste("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/1.Diapause_data/Two_state.nwk",sep=""))
marinebrooted = marineb
mb1 = marineb
mb1$edge.length = c(rep(1,length(mb1$edge.length)))
par(mfrow=c(1,2))
plot(marinebrooted, main="Trait tree from file (1)")
plot(toyTrees$masterTree,main="Manually specified binary tree (1)")
phenvMarine=tree2Paths(marineb, toyTrees)

marineextantforeground = c("Nymphalis_c_album","Nymphalis_polychloros","Nymphalis_urticae","Nymphalis_io","Lysandra_coridon","Plebejus_argus","Aporia_crataegi","Argynnis_bischoffii_washingtonia","Aricia_artaxerxes","Aricia_agestis","Boloria_euphrosyne","Brenthis_hecate","Brenthis_daphne","Brenthis_ino","Coenonympha_glycerion","Cyaniris_semiargus","Erebia_aethiops","Erebia_ligea","Euphydryas_editha","Fabriciana_adippe","Hipparchia_semele","Lasiommata_megera","Lysandra_bellargus","Maniola_hyperantus","Maniola_jurtina","Melanargia_galathea","Melitaea_cinxia","Mellicta_athalia","Oeneis_ivallda","Phengaris_arion","Polyommatus_icarus","Pararge_aegeria","Boloria_selene","Anthocharis_cardamines","Battus_philenor","Celastrina_argiolus","Glaucopsyche_alexis","Leptidea_sinapis","Leptidea_reali","Papilio_xuthus","Papilio_machaon","Pieris_napi","Pieris_mannii","Pieris_rapae","Pieris_brassicae","Hypolimnas_misippus","Pieris_melete")

phenvMarine=foreground2Paths(marineextantforeground, toyTrees, clade="all")
corMarine=correlateWithBinaryPhenotype(mamRERw, phenvMarine, min.sp=10, min.pos=2,
                                       weighted="auto")
head(corMarine[order(corMarine$P),])
corMarine[order(corMarine$P),]
pdf("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/3.HYPHY branch length optimization/5.Plots/period_new_2.pdf", width = 8, height = 5)
plotRers(mamRERw,"period",phenv=phenvMarine)
dev.off()
pdf("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/3.HYPHY branch length optimization/5.Plots/timeless.pdf", width = 8, height = 5)
plotRers(mamRERw,"timeless",phenv=phenvMarine)
dev.off()
pdf("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/3.HYPHY branch length optimization/5.Plots/clock.pdf", width = 8, height = 5)
plotRers(mamRERw,"clock",phenv=phenvMarine)
dev.off()
pdf("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/3.HYPHY branch length optimization/5.Plots/cycle.pdf", width = 8, height = 5)
plotRers(mamRERw,"cycle",phenv=phenvMarine)
dev.off()
pdf("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/3.HYPHY branch length optimization/5.Plots/cry2.pdf", width = 8, height = 5)
plotRers(mamRERw,"cry2",phenv=phenvMarine)
dev.off()
pdf("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/3.HYPHY branch length optimization/5.Plots/cry1a.pdf", width = 8, height = 5)
plotRers(mamRERw,"cry1a",phenv=phenvMarine)
dev.off()
pdf("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/3.HYPHY branch length optimization/5.Plots/cwo.pdf", width = 8, height = 5)
plotRers(mamRERw,"cwo",phenv=phenvMarine)
dev.off()
pdf("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/3.HYPHY branch length optimization/5.Plots/cry1b.pdf", width = 8, height = 5)
plotRers(mamRERw,"cry1b",phenv=phenvMarine)
dev.off()
pdf("G:/My Drive/Circadian Rhythm Genes Project/20.RERconverge/3.HYPHY branch length optimization/5.Plots/timeout.pdf", width = 8, height = 5)
plotRers(mamRERw,"timeout",phenv=phenvMarine)
dev.off()
plotRers(mamRERw,"wg",phenv=phenvMarine)
plotRers(mamRERw,"cad",phenv=phenvMarine)

hist(corMarine$Rho, breaks=15, xlab="Kendall P-value",
     main="P-values for correlation between 200 genes and marine environment")

write.csv(corMarine, 'correlation_2.csv')
