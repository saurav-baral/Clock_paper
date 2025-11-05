# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 11:25:58 2024

@author: sauba
"""
import os
import subprocess

file_location = os.getcwd()
family_group = file_location.split("/")[-1]


location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/15.Timeout_Exon_Analysis/{family_group}"

list_of_species = os.listdir(f"{location}/1.Blast_result")

if "desktop.ini" in list_of_species:
    list_of_species.remove("desktop.ini")

with open(f"{location}/error_in_first_and_last.txt", 'r') as error_file:
    error_list = error_file.readlines()
    
for error in error_list:
    print(error)
    error_split = error.split(",")
    error_species = error_split[0]
    error_exon = error_split[6].split("_")[-1]
    print(error_species, error_exon)

    #making backup
        
    list_of_files_in_species_folder = os.listdir(f"{location}/1.Blast_result/{error_species}")
    
    if "final_coordinates_backup.csv" not in list_of_files_in_species_folder:
        subprocess.run(f'cp "{location}/1.Blast_result/{error_species}/final_coordinates.csv" "{location}/1.Blast_result/{error_species}/final_coordinates_backup.csv"', shell = True)
    
    
    with open(f"{location}/1.Blast_result/{error_species}/final_coordinates_backup.csv", 'r') as coordinate_file:
        coordinate_list = coordinate_file.readlines()
    
    if error_exon == "1":
        
        coordinate_list[1] = error
#         print("".join(coordinate_list))
    else:
        coordinate_list[-1] = error
        
    with open(f"{location}/1.Blast_result/{error_species}/final_coordinates.csv", 'w') as new_coordinate_file:
        new_coordinate_file.write("".join(coordinate_list))