library(RERconverge)

setwd("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/trees_split")
comptrees=readTrees("combined.tree")


mamRERw=RERconverge::getAllResiduals(comptrees, transform = "sqrt", weighted = T, scale = T)

df_t = t(mamRERw)

library(dplyr)
library(cowplot)
library(ggplot2)
library(ggpubr)


diapausing <- c("Anthocharis_cardamines","Aporia_crataegi","Argynnis_bischoffii_washingtonia","Aricia_artaxerxes","Aricia_agestis","Battus_philenor","Boloria_euphrosyne","Boloria_selene","Brenthis_hecate","Brenthis_daphne","Brenthis_ino","Calycopis_cecrops","Celastrina_argiolus","Coenonympha_glycerion","Colias_nastes","Cyaniris_semiargus","Danaus_plexippus","Erebia_aethiops","Erebia_ligea","Euphydryas_editha","Fabriciana_adippe","Glaucopsyche_alexis","Hipparchia_semele","Lasiommata_megera","Leptidea_juvernica","Leptidea_sinapis","Leptidea_reali","Lysandra_coridon","Lysandra_bellargus","Maniola_jurtina","Maniola_hyperantus","Melanargia_galathea","Melitaea_cinxia","Mellicta_athalia","Nymphalis_c-album","Nymphalis_polychloros","Nymphalis_urticae","Nymphalis_io","Oeneis_ivallda","Papilio_xuthus","Papilio_machaon","Papilio_glaucus","Pararge_aegeria","Parnassius_glacialis","Phengaris_arion","Pieris_napi","Pieris_mannii","Pieris_rapae","Pieris_brassicae","Pieris_melete","Plebejus_argus","Polyommatus_icarus","Polyommatus_iphigenia","Colias_palaeno")

non_diapausing <- c("Bicyclus_anynana","Colias_eurytheme","Colias_croceus","Danaus_chrysippus","Dircenna_loreta","Dryadula_phaetusa","Dryas_iulia_moderata","Eueides_isabella","Eurema_hecabe","Heliconius_charithonia","Heliconius_nattereri","Heliconius_sara","Hypolimnas_misippus","Leptophobia_aripa","Mechanitis_messenoides","Mechanitis_mazaeus","Melinaea_menophilus","Melinaea_marsaeus_rileyi","Napeogenes_sylphis","Napeogenes_inachia","Ornithoptera_alexandrae","Ornithoptera_priamus","Papilio_demoleus","Papilio_protenor","Papilio_polytes","Papilio_memnon","Papilio_dardanus_tibullus","Papilio_elwesi","Philaethria_dido","Phoebis_sennae","Teinopalpus_imperialis","Troides_oblongomaculatus","Troides_aeacus","Vanessa_cardui","Vanessa_atalanta","Vanessa_tameamea","Zerene_cesonia")


gene_x <- "cwo.fas.best.fas.ClipKIT"
gene_y <- "period.fas.best.fas.ClipKIT"




plot_df <- as.data.frame(df_t[, c(gene_x, gene_y)])
plot_df <- na.omit(plot_df)



plot_df$species <- rownames(plot_df)


plot_df <- plot_df %>%
  mutate(color_group = case_when(
    species %in% diapausing ~ "Diapausing",
    species %in% non_diapausing ~ "Non-Diapausing",
    TRUE ~ "other"
  ))


plot = ggplot(plot_df, aes_string(x = gene_x, y = gene_y)) +
  geom_point(aes(color = color_group), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_cor(method = "pearson", aes(label = paste0("r = ", ..r..)), size = 5) +
  scale_color_manual(values = c("Diapausing" = "#EF4446", "Non-Diapausing" = "#1A9ACF", "other" = "black")) +
  theme_half_open() +
  background_grid()+
  xlab("CWO") +
  ylab("Period")+
  theme(legend.position="none")


pdf("G:/My Drive/Circadian Rhythm Genes Project/55.Plots_for_Figure_3_coorelations/2.Correlations/Per-cwo.pdf", width = 6, height = 4)
plot
dev.off()




