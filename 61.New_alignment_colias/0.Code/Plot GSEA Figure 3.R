
# setwd("G:/My Drive/Circadian Rhythm Genes Project/53.Microevolutionary_analysis_statistics/2.Pieris_napi/4.Intersections_test/1.For GSEA/2.BUSCO set")

# data = read.table("intersection_set.txt.BP.plotting.csv", sep = ",", header = T)







library(ggplot2)
library(tidyverse)

setwd("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/16.GSEA_Rho-set")
data = read.table("gene_names.txt.BP.plotting.csv", sep = ",", header = T)


sorted_df <-data %>%
  arrange((Enrichment.Score))
sorted_df

sorted_df$Enriched.Ontological.Terms <- factor(sorted_df$Enriched.Ontological.Terms, levels = sorted_df$Enriched.Ontological.Terms)


plot = ggplot(sorted_df, aes(x = Enrichment.Score, y = as.factor(Enriched.Ontological.Terms), size = Number.of.Genes, color = p.value)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.0005), oob = scales::squish) +
  labs(x = "Enrichment Score", y = '', color = "p-value") +
  scale_size_continuous(name = "Gene count")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

pdf("plot_gsea.pdf", width = 8, height = 4.5)
plot
dev.off()


setwd("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/6.GSEA_rho/1.Positive")

data = read.table("gene_names.txt.BP.plotting.csv", sep = ",", header = T)

sorted_df <-data %>%
  arrange((Enrichment.Score))
sorted_df

sorted_df$Enriched.Ontological.Terms <- factor(sorted_df$Enriched.Ontological.Terms, levels = sorted_df$Enriched.Ontological.Terms)

plot = ggplot(sorted_df, aes(x = Enrichment.Score, y = as.factor(Enriched.Ontological.Terms), size = Number.of.Genes, color = p.value)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.0005), oob = scales::squish) +
  labs(x = "Enrichment Score", y = '', color = "p-value") +
  scale_y_discrete(position = "right") +
  scale_x_reverse() +
  scale_size_continuous(name = "Gene count", limits = c(1, 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

pdf("plot_gsea.pdf", width = 8, height = 4.5)
plot
dev.off()







setwd("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/6.GSEA_rho/3.Positive_modi/")

data = read.table("gene_names.txt.BP.plotting.csv", sep = ",", header = T)

sorted_df <-data %>%
  arrange((Enrichment.Score))
sorted_df

sorted_df$Enriched.Ontological.Terms <- factor(sorted_df$Enriched.Ontological.Terms, levels = sorted_df$Enriched.Ontological.Terms)

plot = ggplot(sorted_df, aes(x = Enrichment.Score, y = as.factor(Enriched.Ontological.Terms), size = Number.of.Genes, color = p.value)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.0005), oob = scales::squish) +
  labs(x = "Enrichment Score", y = '', color = "p-value") +
  scale_y_discrete(position = "right") +
  scale_x_reverse() +
  scale_size_continuous(name = "Gene count", limits = c(1, 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

pdf("plot_gsea.pdf", width = 8, height = 4.5)
plot
dev.off()