# -*- coding: utf-8 -*-
"""
Created on Fri May 31 11:35:50 2024

@author: sauba
"""
import os


family_group = "3.Satyrine"
total_exons = 13
location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/10.Cycle_Exon_Analysis/{family_group}"
blast_output_location = f"{location}/1.Blast_result"
species_list = os.listdir(blast_output_location)


if "desktop.ini" in species_list:
    species_list.remove("desktop.ini")

output = "Species"
for i in range(2,total_exons):
    output += f"\tExon_{i}"

for species in species_list:
    output += f"\n{species}"
    query_location = f"{location}/1.Query"
    list_of_query_species = os.listdir(query_location)

    if "desktop.ini" in list_of_query_species:
        list_of_query_species.remove("desktop.ini")

    for query_species in list_of_query_species:
        list_of_processed_exons = os.listdir(f"{blast_output_location}/{species}/Period_gene_genomic_sequence_individual_exon_{query_species}")
#         print(list_of_processed_exons)
        for i in range(2,total_exons):
            if f"Exon_{i}" in list_of_processed_exons:
                output += "\tPresent"
            else:
                output += "\tAbsent"
#         print(species, list_of_query_species)
with open(f"{location}/3.MAFFT_output_check.txt",'w') as out_file:
    out_file.write(output)