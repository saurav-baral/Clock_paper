library(tidyverse)

raw_data = read.table ("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/trees_split/ERC_correlation_full.csv", header = T, sep = ",")
upper_cor_data <- raw_data
upper_cor_data[lower.tri(upper_cor_data,diag = TRUE)] <-NA



cor_data <- upper_cor_data %>%  mutate(colnames=names(raw_data)) %>%  column_to_rownames(var = "colnames") %>% rownames_to_column("row") %>% pivot_longer(cols=-row,names_to = "nodeB", values_to = "ERC_corr",values_drop_na=TRUE)


plot = ggplot(cor_data, aes(x = ERC_corr)) +  geom_histogram(binwidth = 0.05, fill = "grey", color = "black") +   
  labs(title = "ERC correlation histogram", x = "ERC correlations", y = "counts") +     
  geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1) +
  theme_minimal()

plot

pdf("G:/My Drive/Circadian Rhythm Genes Project/55.Plots_for_Figure_3_coorelations/2.Correlations/Historgram_full.pdf", width = 6, height = 4)
plot
dev.off()