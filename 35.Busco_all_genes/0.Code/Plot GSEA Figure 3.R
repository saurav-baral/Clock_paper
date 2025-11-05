
setwd("G:/My Drive/Circadian Rhythm Genes Project/35.Busco_all_genes/11.GSEA_ERC_top_100/GSEA_timeless/")

data = read.table("plot_gsea.csv", sep = ",", header = T)

library(ggplot2)
library(tidyverse)

sorted_df <-data %>%
  arrange((Enrichment.Score))
sorted_df

sorted_df$Enriched.Ontological.Terms <- factor(sorted_df$Enriched.Ontological.Terms, levels = sorted_df$Enriched.Ontological.Terms)

plot = ggplot(sorted_df, aes(x = (Enrichment.Score), y = as.factor(Enriched.Ontological.Terms), size = Number.of.Genes, color = p.value)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = expression("GO terms from analysis of proteins with significant ERC with "*italic("timeless")), x = "Enrichment Score", y = "Enriched Ontological Term")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y=element_text(size=11)) +
  scale_size_continuous(name="Gene count")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("plot_gsea.pdf", width = 10, height = 5)
plot
dev.off()
