library(ggplot2)
library(cowplot)

dataset = read.table("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/7.Manhatten/1.Density_plots_5/pararge_density_plot_x_244_p_0.002.csv", sep = ",", header = T)

xval = 244
p_value = 0.002

p = ggplot(dataset, aes(x = Count)) +
  geom_density(fill = "grey", alpha = 0.5,adjust = 1.5) +
  geom_vline(aes(xintercept = xval),
             color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = xval, y = 0.02, 
           label = paste0("X = ",xval,"\nP-value = ",p_value),
           color = "blue", size = 10, fontface = "bold", vjust = -0.25, angle = -90) +
  
  
  labs(x = "Count", y = "Density")+
  theme_cowplot(12)+
  theme(
    # axis.title.y = element_blank(),   # remove y-axis title
    axis.text.y  = element_blank(),   # remove y-axis labels
    axis.ticks.y = element_blank(),
    # axis.line.y = element_blank(),
    
    axis.title.x = element_text(size = 20),  # increase x-axis label
    # axis.title.y = element_text(size = 20),   # increase y-axis label
    axis.text.x  = element_text(size = 15, face = "bold"), # x-axis tick labels
    # axis.text.y  = element_blank()  # y-axis tick labels
  )

p

pdf("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/7.Manhatten/1.Density_plots_5/pararge_density_plot_x_244_p_0.002.pdf", height = 4, width = 8)
p
dev.off()

