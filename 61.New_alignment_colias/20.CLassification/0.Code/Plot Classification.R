
library(ggplot2)
library(dplyr)
library(cowplot)

# Your species order
species_order <- c(
  'Aricia_agestis', 'Aricia_artaxerxes', 'Polyommatus_icarus', 'Lysandra_bellargus', 'Lysandra_coridon', 'Cyaniris_semiargus', 'Plebejus_argus', 'Celastrina_argiolus', 'Glaucopsyche_alexis', 'Phengaris_arion', 'Calycopis_cecrops', 'Coenonympha_glycerion', 'Hipparchia_semele', 'Oeneis_ivallda', 'Melanargia_galathea', 'Maniola_hyperantus', 'Maniola_jurtina', 'Erebia_aethiops', 'Erebia_ligea', 'Bicyclus_anynana', 'Lasiommata_megera', 'Pararge_aegeria', 'Melitaea_cinxia', 'Mellicta_athalia', 'Euphydryas_editha', 'Hypolimnas_misippus', 'Nymphalis_c-album', 'Nymphalis_polychloros', 'Nymphalis_urticae', 'Nymphalis_io', 'Vanessa_atalanta', 'Vanessa_tameamea', 'Vanessa_cardui', 'Argynnis_bischoffii_washingtonia', 'Fabriciana_adippe', 'Boloria_euphrosyne', 'Boloria_selene', 'Brenthis_daphne', 'Brenthis_ino', 'Brenthis_hecate', 'Dryadula_phaetusa', 'Dryas_iulia_moderata', 'Philaethria_dido', 'Heliconius_charithonia', 'Heliconius_sara', 'Heliconius_nattereri', 'Eueides_isabella', 'Danaus_chrysippus', 'Danaus_plexippus', 'Dircenna_loreta', 'Napeogenes_inachia', 'Napeogenes_sylphis', 'Mechanitis_mazaeus', 'Mechanitis_messenoides', 'Melinaea_marsaeus_rileyi', 'Melinaea_menophilus', 'Anthocharis_cardamines', 'Aporia_crataegi', 'Leptophobia_aripa', 'Pieris_melete', 'Pieris_napi', 'Pieris_mannii', 'Pieris_rapae', 'Pieris_brassicae', 'Colias_eurytheme', 'Colias_palaeno', 'Colias_croceus', 'Zerene_cesonia', 'Phoebis_sennae', 'Eurema_hecabe', 'Leptidea_juvernica', 'Leptidea_reali', 'Leptidea_sinapis', 'Papilio_glaucus', 'Papilio_elwesi', 'Papilio_dardanus_tibullus', 'Papilio_demoleus', 'Papilio_memnon', 'Papilio_polytes', 'Papilio_protenor', 'Papilio_machaon', 'Papilio_xuthus', 'Battus_philenor', 'Ornithoptera_alexandrae', 'Ornithoptera_priamus', 'Troides_aeacus', 'Troides_oblongomaculatus', 'Teinopalpus_imperialis', 'Parnassius_glacialis'
)

df = read.table("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/20.CLassification/test_200.csv", header = T, sep = ",")
# Compute confidence and color

# 
df_long <- df %>%
  pivot_longer(cols = c(True, False, Confused),
               names_to = "Category", values_to = "Proportion")

# Add total test in brackets to species name
df_long <- df_long %>%
  mutate(Species_label = paste0(Species, " (", Total_test, ")"))
  
df_long$Species_label <- factor(df_long$Species_label,
                                levels = rev(paste0(species_order, " (", df$Total_test[match(species_order, df$Species)], ")"))
)


# Plot
plot = ggplot(df_long, aes(x = Proportion, y = Species_label, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("True" = "orange",
                               "False" = "magenta",
                               "Confused" = "lightgreen")) +
  labs(x = "Percentage", y = "Species (Total Tests)",
       fill = "Category") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

# plot = ggplot(df_long, aes(x = Proportion, y = Species_label, fill = Category)) +
#   geom_bar(stat = "identity") +
#   scale_x_continuous(labels = scales::percent) +
#   scale_fill_manual(values = c("True" = "red",
#                                "False" = "blue",
#                                "Confused" = "lightgreen")) +
#   labs(x = "Percentage", y = "Species (Total Tests)",
#        fill = "Category") +
#   theme_minimal() +
#   theme(axis.text.y = element_text(size = 10))

pdf("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/20.CLassification/test_200.pdf", height = 15, width = 8)
plot
dev.off()
