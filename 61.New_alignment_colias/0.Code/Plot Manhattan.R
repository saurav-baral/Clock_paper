library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)


# Variables ====
mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette
mypalette <- colorRampPalette(c("grey70", "black"))(22)

mypalette <- c("grey70", "black")
mysnps <- c("rs11801961","rs116558464","rs61703161") # snps to highlight
sig = 0.1 # significant threshold line
sugg = 0.2 # suggestive threshold line

gg.manhattan <- function(df, threshold, hlight, col, ylims, title){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    #mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate(is_highlight = ifelse(!is.na(SNP) & SNP != "", "yes", "no")) %>%
    #mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
    mutate(is_annotate=ifelse(!is.na(Common) & Common != "", "yes", "no"))
  
  # get chromosome center positions for x-axis
  # axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  # axisdf <- df.tmp %>% summarize(center = (max(BPcum) + min(BPcum)) / 2, .by = CHR)
  # axisdf <- df.tmp %>%
  #   group_by(CHR) %>%
  #   summarize(center = (max(BPcum) + min(BPcum)) / 2, .groups = "drop")
  
  df.tmp$CHR <- as.character(df.tmp$CHR)
  
  axisdf <- df.tmp %>%
    group_by(CHR) %>%
    summarise(center = (max(BPcum) + min(BPcum)) / 2)
  
  #ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
  ggplot(df.tmp, aes(x=BPcum, y=(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1) +
    geom_point(data = subset(df.tmp, is_highlight == "yes"), color = "red", size = 1) +
    # scale_color_manual(values = rep(col, 22 )) +
    scale_color_manual(values = rep(col, length(unique(dataset$CHR))))+
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    labs(y = "Fst") +
    
    # # add genome-wide sig and sugg lines
    # geom_hline(yintercept = -log10(sig)) +
    # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    # geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    geom_label_repel(
      data = df.tmp %>%
        filter(is_annotate == "yes") %>%
        distinct(SNP, .keep_all = TRUE),
      aes(label = as.factor(SNP)),
      size = 1.5,
      force = 20,
      # nudge_y = -0.05,
      max.overlaps = Inf,
      alpha = 0.7
      
    ) +
    
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme_cowplot(12) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, vjust = .5, hjust =1)
    )
    
}





dataset = read.table("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/7.Manhatten/dataset_5_lasiommata.csv", header = T, sep = ",")


# Extract unique chromosome names

chr_levels <- unique(dataset$CHR)

# Separate numeric and non-numeric chromosome names
numeric_chr <- suppressWarnings(as.numeric(as.character(chr_levels)))
is_numeric <- !is.na(numeric_chr)

# Sort numeric chromosomes and the rest separately
sorted_chr <- c(
  as.character(sort(numeric_chr[is_numeric])),
  sort(as.character(chr_levels[!is_numeric]))
)

# Apply sorted factor levels
dataset$CHR <- factor(dataset$CHR, levels = sorted_chr)


pdf("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/7.Manhatten/dataset_5_lasiommata.pdf", width = 11, height = 2.5)
gg.manhattan(dataset, threshold=2, hlight=character(0),col=mypalette, ylims=c(0.0,1.1), title="")
dev.off()

#print(gg.manhattan(loco_asm_fev, sugg, NA, c(1.5,9), reds.c, "Pre-treatment FEV1 % Predicted - MLMA (LOCO) GWAS Results\n"))
