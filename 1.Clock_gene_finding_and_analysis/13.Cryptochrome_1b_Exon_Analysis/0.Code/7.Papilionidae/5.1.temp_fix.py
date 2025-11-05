# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 12:24:45 2024

@author: sauba
"""
import os
import subprocess

file_location = os.getcwd()
family_group = file_location.split("/")[-1]


annotated_genome_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/13.Cryptochrome_1b_Exon_Analysis/{family_group}/1.Blast_result"

species_list = os.listdir(annotated_genome_location)

if "desktop.ini" in species_list:
    species_list.remove("desktop.ini")
    
for species in species_list:
    subprocess.run(f'cp "{annotated_genome_location}/{species}/{species}_coordinates_old.csv" "{annotated_genome_location}/{species}/final_coordinates.csv"', shell=True)