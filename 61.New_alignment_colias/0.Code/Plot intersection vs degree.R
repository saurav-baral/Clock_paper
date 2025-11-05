library(ggplot2)

data = read.table("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/14.Plot degree vs intersction/intersection_countvsdegree.csv", sep = ",", header = T)

pdf("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/14.Plot degree vs intersction/intersection_countvsdegree.pdf", width = 8, height = 4.5)

ggplot(data, aes(x = as.factor(Intersection), y = Degree)) +
  geom_boxplot()+
  geom_jitter()

dev.off()
