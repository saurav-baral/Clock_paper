# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 11:04:39 2024

@author: sauba
"""
import os
import re
import io

list_of_species = ["Hestina_assimilis","Charaxes_varanes","Danaus_chrysippus","Danaus_plexippus","Dircenna_loreta","Mechanitis_mazaeus","Mechanitis_messenoides","Melinaea_marsaeus_rileyi","Melinaea_menophilus","Napeogenes_inachia","Napeogenes_sylphis","Argynnis_bischoffii_washingtonia","Boloria_euphrosyne","Boloria_selene","Brenthis_daphne","Brenthis_hecate","Brenthis_ino","Dryadula_phaetusa","Dryas_iulia_moderata","Eueides_isabella","Fabriciana_adippe","Heliconius_charithonia","Heliconius_nattereri","Heliconius_sara","Philaethria_dido","Limenitis_camilla","Euphydryas_editha","Hypolimnas_misippus","Melitaea_cinxia","Mellicta_athalia","Nymphalis_c-album","Nymphalis_io","Nymphalis_polychloros","Nymphalis_urticae","Vanessa_atalanta","Vanessa_cardui","Vanessa_tameamea"]
query_species = "Danaus_plexippus"

file_location = os.getcwd()
family_group = file_location.split("/")[-1]

def get_genome_file(genome_location,species):
    list_of_files_in_genome_folder = os.listdir(f"{genome_location}/{species}")
    for file in list_of_files_in_genome_folder:
        if file.endswith("_genomic.fna"):
            genome_file = file
    return(genome_file)   

gene_folder = "11.CWO_Exon_Analysis"


query_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/11.CWO_Exon_Analysis/{family_group}/1.Query"
query_transcript_list = os.listdir(f"{query_location}/{query_species}")
if "desktop.ini" in query_transcript_list:
    query_transcript_list.remove("desktop.ini")

query_transcript = query_transcript_list[0]

list_of_genomes = os.listdir("/mnt/griffin/saubar/Genomes_2023-12-26")
query_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/11.CWO_Exon_Analysis/{family_group}/1.Query/{query_species}/{query_transcript}/query_protein.fa"

genome_file = ''
output = ''
for species in list_of_species:
    species = species.replace(" ", "_")
    # if "
    # species = species_folder.split(".")[1]
    blast_output_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/11.CWO_Exon_Analysis/{family_group}/1.Blast_result"
    try:
        os.mkdir(blast_output_location)
    except:
        pass
    
    blast_output_folders = os.listdir(blast_output_location)
    if species not in blast_output_folders:
        os.mkdir(f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/11.CWO_Exon_Analysis/{family_group}/1.Blast_result/{species}")
    
    genome_location = "/mnt/griffin/saubar/Genomes_2023-12-26"
#     for genome in list_of_genomes:
#         if species == genome:
#             files_in_genome_folder = os.listdir(f"/mnt/f/Genomes_2023-12-26/{genome}")
#             for file in files_in_genome_folder:
#                 if re.search('_genomic.fna$', file):
#                     genome_file = file
#                     break
    genome_file = get_genome_file(genome_location,species)
    
    output += f'echo {species}\ntblastn -seg no -query "{query_location}" -db "{genome_location}/{species}/{genome_file}" -num_alignments 3 -out "{blast_output_location}/{species}/{species}_blast_out.htm" -html \n\n tblastn -seg no -query "{query_location}" -db "{genome_location}/{species}/{genome_file}" -num_alignments 3 -out "{blast_output_location}/{species}/{species}_blast_out.txt"\n\n'
    
with io.open(f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/11.CWO_Exon_Analysis/{family_group}/tblastn_{query_species}.sh", "w", newline='\n') as outfile:
    outfile.write(output)


