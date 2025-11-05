##Pipeline to compare the covariation of physically interacting domains versus non-interacting domains
##Updated by Jordan Little 9/28/2023

##package dependencies
library(PRROC)
library(phangorn)
library(plyr)
library(reshape2)
library(RERconverge)


prep_complexes = function(name, full_mat,complex_association, tp_list){
  ##make complex matrix symmetrical and only take the bottom half so we don't get
  ##duplicate edges
  complex_list = as.character(complex_association$domain[complex_association$complex == name])
  tp_list = tp_list[tp_list$GENEA %in% complex_list,]
  
  complex_mat = full_mat[complex_list, complex_list]
  complex = make_symmetric(complex_mat)
  complex[upper.tri(complex)] = NA
  
  raw_genes = unique(sapply(strsplit(complex,'_'),'[',1))
  ##set self x self cells to NA
  for(i in 1:length(raw_genes)){
    complex[grepl(paste(raw_genes[i],"_"), rownames(complex)), grepl(paste(raw_genes[i],"_"), colnames(complex))] = NA
  }
  
  
  domain_edges = reshape2::melt(complex, na.rm = TRUE)
  colnames(domain_edges) = c("GENEA","GENEB","ftERC")
  
  
  ##Turn the true positive list into a matrix so that we can smoosh it with the 
  ##erc edges 
  tp_mat = matrix(data = NA, nrow = length(complex_list), ncol = length(complex_list))
  rownames(tp_mat) = complex_list
  colnames(tp_mat) = complex_list
  
  for(i in 1:nrow(tp_list)){
    genea = tp_list[,1][i]
    geneb = tp_list[,2][i]
    tp_mat[genea, geneb] = 1
  }
  
  ##make the true positive matrix symmetrical and only take the bottom half so we
  ##don't have duplicate edges
  tp_mat[is.na(tp_mat)] = 0
  tp_mat[is.na(domain_mat)] = NA
  tp_mat = make_symmetric(tp_mat)
  tp_mat[upper.tri(tp_mat)] = NA
  tp_edges = reshape2::melt(tp_mat, na.rm = TRUE)
  colnames(tp_edges) = c("GENEA","GENEB","true_positive")
  
  ##merge the two edgelists
  domain_tp = merge(domain_edges, tp_edges)
  
  domain_tp
}

plot_roc_domains_only = function(edgelist, title){
  res=edgelist[order(edgelist$ftERC,decreasing = T),]
  
  roc_d=roc.curve(scores.class0 =res[res[,4]=="1",]$ERC,
                  scores.class1 =res[res[,4]=="0",]$ERC,
                  curve=T,)
  plot(roc_d, auc.main = TRUE,main = paste(title, "domain interactions"), color = 'black')
  roc_d
}

plot_all_complexes = function(ROC){
  plot.new()
  par(bg = "#a9a9a9")
  
  par(mar = c(5,4,4,11), xpd = TRUE) ##these dimensions have to be changed based on the plot window size in Rstudio
  all_roc = plot(ROC[[1]], main = "Complex ROC curves",xlab = "False positive rate", max.plot =TRUE, min.plot = TRUE, rand.plot = TRUE, fill.area = T, color =  colors[9], auc.main = FALSE)
  aucs = as.matrix(c(sapply(ROC,'[[',2)))
  aucs = format(round(aucs, 3))
  names = rownames(aucs)
  list_for_legend = mapply(function(x,y) paste(x,"(", y,")", sep = ""), names, aucs) 
  legend("topright", inset = c(-0.85,0.1), title = "Complex (ROC-AUC)",legend = list_for_legend, col = colors, cex=0.8, lty = 1, lwd = 3)
  for(i in 2:length(ROC)){
    plot(ROC[[i]], color = colors[i], add = TRUE)
    op = par(cex = 1.7)
  }

  
}

##get the permutation p-value and False positive rate for each complex
pairwise_rank_permutation = function(erc_mat, complex_edgelist, perms = 1000){
  complex_tp = complex_edgelist[,c("GENEA","GENEB")][complex_edgelist$true_positive ==1]
  complex_domains = unique(as.character(c(complex_edgelist$GENEA, complex_edgelist$GENEB)))
  complex_mat = make_symmetric(erc_mat[complex_domains, complex_domains])
  complex_pval = list()
  obs = NA
  for(true_pos in 1:nrow(complex_tp)){
    genea = sapply(strsplit(as.character(complex_tp$GENEA[true_pos]),'_'),'[',1)
    geneb = sapply(strsplit(as.character(complex_tp$GENEB[true_pos]),'_'),'[',1)
    if(length(complex_pval[[paste(genea,geneb)]]) == 3){
      obs = obs
      cnt = cnt + 1
    }else{
      obs = NA
      cnt = 1
    }
    rows = rownames(complex_mat)[grep(paste(genea,"_",sep=''), rownames(complex_mat))]
    cols = rownames(complex_mat)[grep(paste(geneb,"_",sep=''), rownames(complex_mat))]
    pair_mat = as.matrix(complex_mat[rows,cols])
    
    if(sum(rownames(pair_mat) == rows) != length(rows)){
      pair_mat = t(pair_mat)
      rownames(pair_mat) = rows
      colnames(pair_mat) = cols
    }

    pair_edgelist = reshape2::melt(pair_mat, na.rm = TRUE)
    if(nrow(pair_edgelist) <1){
      next
    }
    pair_edgelist$true_pos = 0
    pair_edgelist = pair_edgelist[order(pair_edgelist$value, decreasing = TRUE),]
    row.names(pair_edgelist) = NULL
    tp = which(pair_edgelist$Var1 == as.character(complex_tp$GENEA[true_pos]) & pair_edgelist$Var2 == as.character(complex_tp$GENEB[true_pos]))
    if(length(tp) ==0){
      tp = which(pair_edgelist$Var1 == as.character(complex_tp$GENEB[true_pos]) & pair_edgelist$Var2 == as.character(complex_tp$GENEA[true_pos]))
    }
    
    pair_edgelist$true_pos[tp] = 1
    temp = (1-tp)/(1-nrow(pair_edgelist))
    if(length(temp) == 0){
      next
    }
    null = c()
    m = c()
    while (length(null) < perms) {
      s <- sample(x=1:nrow(pair_edgelist) , cnt)
      for(k in 1:cnt){
        m[k] = (1-s[k])/(1-nrow(pair_edgelist))
      }
      m = mean(m)
      if (is.nan(m)) next
      null <- c(null , m)
    }
    
    obs = mean(c(temp, obs), na.rm = TRUE)
    pval <- sum(as.numeric(obs >= null)) / perms
    print(pval)
    complex_pval[[paste(genea,geneb)]] = list( obs=obs , p=pval , null=null)
  }
  complex_pval
}

complex_fpr_permutations = function(complex_fpr, perms = 1000){
  obs = mean(sapply(complex_fpr, "[[",1))
  null = rowMeans(as.matrix(sapply(complex_fpr, "[[",3)))
  pval = sum(as.numeric(obs>=null)) / perms
  complex_pval = list(obs=obs, p = pval, null=null)
  complex_pval
}

setwd("C:/Users/Saurav Baral/Desktop/Busco_related/20.RERConverge/")
estimatePhangornTreeAll(alndir = "1.Alignment", treefile = "combined_file.fas.treefile_no_branch_length.nwk", output.file = "trees_with_branch_length.tre")
