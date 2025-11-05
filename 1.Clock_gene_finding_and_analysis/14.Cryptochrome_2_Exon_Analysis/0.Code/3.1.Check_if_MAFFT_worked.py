# -*- coding: utf-8 -*-
"""
Created on Fri May 31 11:35:50 2024

@author: sauba
"""
import os

list_of_family_groups = ["3.Satyrine","4.Pierinae","5.Coliadinae","6.Heliconiinae_Danainae_Nymphalinae","7.Papilionidae","8.Lycaenidae"]
output = ''
gene_folder = "14.Cryptochrome_2_Exon_Analysis"
total_exons = 9
output += "Family\tSpecies\tQuery"
for i in range(2,total_exons):
    output += f"\tExon_{i}"
    
for family_group in list_of_family_groups:
#family_group = "3.Satyrine"
    
    location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/{gene_folder}/{family_group}"
    blast_output_location = f"{location}/1.Blast_result"
    species_list = os.listdir(blast_output_location)


    if "desktop.ini" in species_list:
        species_list.remove("desktop.ini")

    

    for species in species_list:
        
        query_location = f"{location}/1.Query"
        list_of_query_species = os.listdir(query_location)

        if "desktop.ini" in list_of_query_species:
            list_of_query_species.remove("desktop.ini")

        for query_species in list_of_query_species:
            output += f"\n{family_group}\t{species}\t{query_species}"
            try:
                list_of_processed_exons = os.listdir(f"{blast_output_location}/{species}/Period_gene_genomic_sequence_individual_exon_{query_species}")
            except:
                list_of_processed_exons = []
    #         print(list_of_processed_exons)
            for i in range(2,total_exons):
                if f"Exon_{i}" in list_of_processed_exons:
                    output += "\tPresent"
                else:
                    output += "\tAbsent"
    #         print(species, list_of_query_species)
print(output)
with open(f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/{gene_folder}/3.MAFFT_output_check.txt",'w') as out_file:
    out_file.write(output)