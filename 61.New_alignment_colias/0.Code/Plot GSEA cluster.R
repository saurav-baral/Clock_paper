library(ggplot2)
library(tidyverse)
folder_list <- list.dirs(path = "C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/17.GSEA_Rho-set_clusters/", full.names = TRUE, recursive = FALSE)

for (folder in folder_list) {
  # setwd ("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/17.GSEA_Rho-set_clusters/Cluster_0")
  setwd(folder)
  print(folder)
  cluster_number <- sub("Cluster_", "", basename(folder))
  # print(cluster_number)
  data = read.table("gene_names.txt.BP.plotting.csv", sep = ",", header = T)
  
  
  sorted_df <-data %>%
    arrange((Enrichment.Score))
  sorted_df
  
  sorted_df$Enriched.Ontological.Terms <- factor(sorted_df$Enriched.Ontological.Terms, levels = sorted_df$Enriched.Ontological.Terms)
  
  
  plot = ggplot(sorted_df, aes(x = Enrichment.Score, y = as.factor(Enriched.Ontological.Terms), size = Number.of.Genes, color = p.value)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.0005), oob = scales::squish) +
    labs(x = paste("Enrichment Score Cluster", cluster_number), y = '', color = "p-value") +
    scale_size_continuous(name = "Gene count")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

  pdf("plot_gsea.pdf", width = 8, height = 4.5)
  print(plot)
  dev.off()
}
