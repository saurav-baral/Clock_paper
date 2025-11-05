library(ggplot2)
library(cowplot)
omega_values = read.table("G:/My Drive/Circadian Rhythm Genes Project/21.Omega_value_results/omega_values.txt", header = T, sep = "\t")
omega_values$Gene <- factor(omega_values$Gene, levels = c("Period","Timeless","CLOCK/BMAL1", "Cycle","CWO","Cry1a","Cry1b","Cry2","Timeout"))
plot = ggplot(data = omega_values) +
  geom_boxplot(aes(x = Gene, y = Omega.dN.dS.))+
  geom_jitter(aes(x = Gene, y = Omega.dN.dS.,col = Family), width = 0.2,alpha = 0.8)+
  ylim(0, 1)+
  xlab("Clock Genes")+
  ylab("Omega (dN/dS)")+ 
  theme_cowplot(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
plot
pdf("G:/My Drive/Circadian Rhythm Genes Project/21.Omega_value_results/combined_without_wg.pdf", width = 8, height = 5)
plot
dev.off()
plot = ggplot(data = omega_values) +
  geom_boxplot(aes(x = Gene, y = Omega.dN.dS.,col = Family))+
  
  ylim(0, 1)+
  xlab("Clock Genes")+
  ylab("Omega (dN/dS)")+
  theme_cowplot(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
plot
pdf("G:/My Drive/Circadian Rhythm Genes Project/21.Omega_value_results/family_specific_without_wg.pdf", width = 8, height = 5)
plot
dev.off()
