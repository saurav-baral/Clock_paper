dataset = read.table("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/22.Upset_plot_deconstructed_for_degree/group_degree.csv", header= T, sep = ',')

library(ggplot2)
library(cowplot)

plot = ggplot (data = dataset) +
  geom_boxplot(aes(x = group, y = degree))+
  geom_jitter(aes(x = group, y = degree, col =group ))+
  theme_cowplot(12) + background_grid()

pdf("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/22.Upset_plot_deconstructed_for_degree/group_degree.pdf", height = 6, width = 10)
plot
dev.off()
