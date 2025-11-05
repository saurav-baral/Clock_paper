dataset = read.table("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/22.Upset_plot_deconstructed_for_degree/3.Intersection_genes_upset_cluster/Cluster.csv", header= T, sep = ',')

library(ggplot2)
library(cowplot)

dataset$Upset.class <- factor(dataset$Upset.class,
                              levels = c("All three", "L. megera & P. aegeria", "L. megera & P. napi","P. napi & P. aegeria","Lasiommata megera", "Pararge aegeria", "Pieris napi"))

plot = ggplot (data = dataset) +
  geom_bar(aes(x = as.factor(Cluster), fill = Upset.class),position = "dodge")+
  xlab("Clusters")+
  ylab("Gene count")+
  
  theme_cowplot(12) + background_grid()

pdf("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/22.Upset_plot_deconstructed_for_degree/Cluster.pdf", height = 6, width = 15)
plot
dev.off()

