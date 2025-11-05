library(ggplot2)
data = read.table("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/18.Cluster_in_each_species/intersection_cluster.csv", sep = ',', header = T)

pdf("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/18.Cluster_in_each_species/plot_count.pdf", width = 8, height = 4.5)

ggplot(data,aes(x = as.factor(Cluster), y = count, fill = Species))+
  geom_col(position = "dodge")+
  theme_minimal() +
  labs(title = "Counts per Cluster per Group",
       x = "Cluster", y = "Gene Count")

dev.off()