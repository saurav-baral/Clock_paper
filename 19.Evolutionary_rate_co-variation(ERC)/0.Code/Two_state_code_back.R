setwd("G:/My Drive/Circadian Rhythm Genes Project/19.Evolutionary_rate_co-variation(ERC)/2.Two_state")

comparison_table = read.table("ERC_correlation_diapause_nd_difference_0.5.csv", header = T, sep = ",")


library(ggrepel)
library(tidyverse)

threshold=0.2398
# Count the number of values greater than the threshold
count_above_threshold <- sum(data$Rho > threshold)
counts <- length(data$Rho)

length(unique(comparison_table$Gene2))


plot = ggplot(comparison_table, aes(x = Dia.NDia)) +
  geom_histogram(binwidth = 0.005, fill = "skyblue", color = "black")   # Histogram
  geom_vline(xintercept = 0.4, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = 0.7,y = Inf, label = paste("> 0.4"), hjust = 1.1, vjust = 1.1, color = "blue", size = 5)+
  labs(title = "Histogram of ERC values", x = "Values", y = "Frequency") +
  theme_minimal()  # Optional, sets a minimal theme
plot

p3 <- comparison_table %>% 
  ggplot(aes(x=Dia.NDia,y=abs(Dia.NDia))) + 
  # xlim(0, 1000) +  # Set x limits from 0 to 6
  # ylim(0, 1000) +  # Set y limits from 0 to 6
  geom_point()   # Color based on condition
  
p3

scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + 
  geom_text_repel(
    data = subset(data, P < 0.005),  # Only label points with p < 0.05
    aes(label = genes),  # Label with gene names
    vjust = -0.5,  # Position text slightly above points
    size = 3,  # Adjust text size
    color = "black"  # Set label color
  ) +
  labs(title = "Scatter Plot with Gene Labels", x = "Rho", y = "-log10(P-value)") +
  geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1) +
  # geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +  # Add threshold line
  annotate("text", x = Inf, y = Inf, label = paste("Count >", threshold, "=", count_above_threshold,";",count_above_threshold,"/total=",round(count_above_threshold/counts,4)), 
           hjust = 1.1, vjust = 1.1, color = "blue", size = 5)  + # Annotate count as text
  theme_minimal()

ggMarginal(p3, type = "histogram")