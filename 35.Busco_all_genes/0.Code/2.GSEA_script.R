#####

# This is a  script to take in two files that have been pre processed
# 1) you need an annotation file that has a list of all genes tab then 
# comma separated GO terms
# 2) you need a list of the candidate genes from that list you want to assess 
# for being enriched compared to the genome.
#
# NOTE: these two files must have the following headers added to them.
# annotation=japonicus_aa.fasta.emapper.annotations
# candidate_file=japponicus_run6above95
# cut -f1,6 $annotation > $annotation.tsv
# echo 'Parent  GO_term' | cat - "$annotation".tsv > "$annotation".header.tsv
# echo 'geneid' | cat - "$candidate_file" > "$candidate_file".header.tsv
#
# the resulting *header.tsv files are then used in the script below
# 1st argument is the annotation file 
# 2nd argument is the candidate set list 
#
# note that the current setting for GO node size is 5, which 
# is hard coded below for more robust results.
#
# example run
# Rscript GSEA_run_script.R gifuensis_aa.fasta.emapper.annotations.header.tsv gifuensis_run6above95.header.tsv

# descriptions
# http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html

# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

library("topGO")
library("Rgraphviz")
library(openxlsx)

##########
# set variables
# making general script 
# Rscript myscript.R batch.csv
# and invoke these in the myscript.R
args <- commandArgs(TRUE)
# dataset <- read.table(args[1],header=FALSE,sep=",",skip=1)

annotations=args[1]
candidate_list=args[2]
# annotations="MM_90qfcby5.emapper.annotations.header.tsv"
# candidate_list="candidate_genes_table.header.tsv"
# here I will be only analyzing GO terms with at least 5 members,
# as this yield more stable results.
node_size=5

#### CORE BP #####
GO_category="BP"
geneID2GO <- readMappings(file = annotations)  
geneUniverse <- names(geneID2GO) 

genesOfInterest.bv <- read.table(candidate_list,header=TRUE)

genesOfInterest.bv <- as.character(genesOfInterest.bv$geneid) 
geneList.bv <- factor(as.integer(geneUniverse %in% genesOfInterest.bv))
names(geneList.bv) <- geneUniverse

myGOdata.bv <- new("topGOdata", description="Candidate genes", ontology=GO_category, allGenes=geneList.bv,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = node_size)

# STATS#
# each GO term is tested independently, not taking the GO hierarchy into account
resultClassic <- runTest(myGOdata.bv, algorithm="classic", statistic="fisher")
# elim method processes the GO terms by traversing the GO hierarchy from bottom to top, 
# ie. it first assesses the most specific (bottom-most) GO terms, and proceeds later 
# to more general (higher) GO terms. When it assesses a higher (more general) GO term, 
# it discards any genes that are annotated with significantly enriched descendant 
# GO terms (considered significant using a pre-defined P-value threshold). 
# This method does tend to miss some true positives at higher (more general) 
# levels of the GO hierarchy.
resultElim <- runTest(myGOdata.bv, algorithm="elim", statistic="fisher")
# weight01 this is the default method used by TopGO, and is a mixture of the 'elim' and 'weight' methods
resultTopgo <- runTest(myGOdata.bv, algorithm="weight01", statistic="fisher")
# when assessing a GO term, it takes into accoount the annotation of terms to the current term's parents, 
# and so reduces false positives due to the inheritance problem
resultParentchild <- runTest(myGOdata.bv, algorithm="parentchild", statistic="fisher")


# see how many results we get where weight01 gives a P-value <= 0.001:
mysummary <- summary(attributes(resultTopgo)$score <= 0.1)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001

allRes <- GenTable(myGOdata.bv, 
				   classicFisher = resultClassic, 
				   elimFisher = resultElim, 
				   topgoFisher = resultTopgo, 
				   parentchildFisher = resultParentchild, 
				   orderBy = "parentchildFisher", ranksOf = "classicFisher", topNodes = numsignif)


# write output
printGraph(myGOdata.bv, resultClassic, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultClassic", sep=""), useInfo = "all", pdfSW = TRUE)
printGraph(myGOdata.bv, resultTopgo, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultTopGo", sep=""), useInfo = "all", pdfSW = TRUE)
printGraph(myGOdata.bv, resultParentchild, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultParentchild", sep=""), useInfo = "all", pdfSW = TRUE)

write.table(allRes[,c(1,8)], file=paste(candidate_list, ".",GO_category,".GSEA_result_elimFisher.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(allRes[,c(1,7)], file=paste(candidate_list, ".",GO_category,".GSEA_result_classicFisher.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(allRes[,c(1,10)], file=paste(candidate_list, ".",GO_category,".GSEA_result_parentchild.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(allRes, file=paste(candidate_list, ".",GO_category,".GSEA_result.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.xlsx(allRes, file = paste(candidate_list, ".",GO_category,".GSEA_result.xlsx", sep=""), borders = "rows")

#### CORE MF #####
# new general run
objects() # clear all objects definded the previous run

##########
# set variables
# making general script 
# Rscript myscript.R batch.csv
# and invoke these in the myscript.R
rm(myGOdata.bv,allRes,resultClassic,resultElim,resultTopgo,resultParentchild,GO_category) #clears all variables
args <- commandArgs(TRUE)
# dataset <- read.table(args[1],header=FALSE,sep=",",skip=1)

annotations=args[1]
candidate_list=args[2]
# here I will be only analyzing GO terms with at least 5 members,
# as this yield more stable results.
node_size=5

#### CORE MF #####
GO_category="MF"
geneID2GO <- readMappings(file = annotations)  
geneUniverse <- names(geneID2GO) 

genesOfInterest.bv <- read.table(candidate_list,header=TRUE)

genesOfInterest.bv <- as.character(genesOfInterest.bv$geneid) 
geneList.bv <- factor(as.integer(geneUniverse %in% genesOfInterest.bv))
names(geneList.bv) <- geneUniverse

myGOdata.bv <- new("topGOdata", description="Candidate genes", ontology=GO_category, allGenes=geneList.bv,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = node_size)

# STATS#
# each GO term is tested independently, not taking the GO hierarchy into account
resultClassic <- runTest(myGOdata.bv, algorithm="classic", statistic="fisher")
# elim method processes the GO terms by traversing the GO hierarchy from bottom to top, 
# ie. it first assesses the most specific (bottom-most) GO terms, and proceeds later 
# to more general (higher) GO terms. When it assesses a higher (more general) GO term, 
# it discards any genes that are annotated with significantly enriched descendant 
# GO terms (considered significant using a pre-defined P-value threshold). 
# This method does tend to miss some true positives at higher (more general) 
# levels of the GO hierarchy.
resultElim <- runTest(myGOdata.bv, algorithm="elim", statistic="fisher")
# weight01 this is the default method used by TopGO, and is a mixture of the 'elim' and 'weight' methods
resultTopgo <- runTest(myGOdata.bv, algorithm="weight01", statistic="fisher")
# when assessing a GO term, it takes into accoount the annotation of terms to the current term's parents, 
# and so reduces false positives due to the inheritance problem
resultParentchild <- runTest(myGOdata.bv, algorithm="parentchild", statistic="fisher")


# see how many results we get where weight01 gives a P-value <= 0.001:
mysummary <- summary(attributes(resultTopgo)$score <= 0.1)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001

allRes <- GenTable(myGOdata.bv, 
				   classicFisher = resultClassic, 
				   elimFisher = resultElim, 
				   topgoFisher = resultTopgo, 
				   parentchildFisher = resultParentchild, 
				   orderBy = "parentchildFisher", ranksOf = "classicFisher", topNodes = numsignif)


# write output
printGraph(myGOdata.bv, resultClassic, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultClassic", sep=""), useInfo = "all", pdfSW = TRUE)
printGraph(myGOdata.bv, resultTopgo, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultTopGo", sep=""), useInfo = "all", pdfSW = TRUE)
printGraph(myGOdata.bv, resultParentchild, firstSigNodes = 5, fn.prefix = paste(candidate_list, ".",GO_category,".GSEA_graph_resultParentchild", sep=""), useInfo = "all", pdfSW = TRUE)

write.table(allRes[,c(1,8)], file=paste(candidate_list, ".",GO_category,".GSEA_result.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(allRes[,c(1,7)], file=paste(candidate_list, ".",GO_category,".GSEA_result_classicFisher.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(allRes[,c(1,10)], file=paste(candidate_list, ".",GO_category,".GSEA_result_parentchild.REVIGO.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(allRes, file=paste(candidate_list, ".",GO_category,".GSEA_result.tsv", sep=""), sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.xlsx(allRes, file = paste(candidate_list, ".",GO_category,".GSEA_result.xlsx", sep=""), borders = "rows")
