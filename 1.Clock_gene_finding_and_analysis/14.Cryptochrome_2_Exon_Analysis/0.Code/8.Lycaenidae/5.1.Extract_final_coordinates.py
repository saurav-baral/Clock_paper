# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 10:11:20 2024

@author: sauba
"""
import os

file_location = os.getcwd()
family_group = file_location.split("/")[-1]




annotated_genome_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/14.Cryptochrome_2_Exon_Analysis/{family_group}/1.Blast_result"

species_list = os.listdir(annotated_genome_location)

if "desktop.ini" in species_list:
    species_list.remove("desktop.ini")

for species in species_list:
    print(species)
    blast_location = annotated_genome_location
    
    list_of_files = os.listdir(f"{blast_location}/{species}")
    query_species_list = []
    coordinate_detail_dictionary = {}
    coordinate_info_dictionary = {}
    scaffold = ''
    for file in list_of_files:
        if file.endswith(".csv") and len(file.split("_")) >= 5 and not(file.endswith("_old.csv")) and not(file.endswith("_old_backup.csv")):
#             print(file)
            if "." in file[:-4]:
                query_species = f"{file.split('_')[-2].split('.')[1]}_{file.split('_')[-1].split('.')[0]}"
            else:
                query_species = f"{file.split('_')[-2]}_{file.split('_')[-1].split('.')[0]}"
            query_species_list.append(query_species)
            coordinate_detail_dictionary[query_species] = {}
            with open(f"{blast_location}/{species}/{file}", 'r') as coordinate_file_open:
                coordinate_file_lines_list  =  coordinate_file_open.readlines()
            if scaffold == '':
                scaffold = coordinate_file_lines_list[1].split(",")[1]
    
            for i in range(1,len(coordinate_file_lines_list)):
                
                coordinate_file_lines_split = coordinate_file_lines_list[i].split(",")
                start,stop = coordinate_file_lines_split[2],coordinate_file_lines_split[3]
                exon = f"Exon_{i}"
                if exon not in coordinate_info_dictionary:
                    coordinate_info_dictionary[exon] = {}
                coordinate_info_dictionary[exon][query_species] = [start,stop]
                coordinate_detail_dictionary[query_species][exon] = coordinate_file_lines_list[i]
            
            # break
            # print(query_species)
            # print(file)
#     print(coordinate_info_dictionary)
#     print(coordinate_detail_dictionary)
    # assert False
    output = coordinate_file_lines_list[0]
    for exon,species_specific_coordinates in coordinate_info_dictionary.items():
        # print(exon,species_specific_coordinates)
        coordinates_list = []
        previous_coordinates = ""
        match = 1
        for query_species, coordinates in coordinate_info_dictionary[exon].items():
            coordinates_list.append([query_species,coordinates])
            if previous_coordinates != '':
                if coordinates != previous_coordinates:
                    match = 0
                    
            previous_coordinates = coordinates
       
        
        
        if match == 0:
            print(exon,"mismatch")
            print(coordinates_list)
                
            while True:
                choice = input("Proceed?")
                if choice[0].lower() == "y" or choice[0].lower() == "n":
                    break
            
            while True:
                if choice[0].lower() == "y":
                    alignment_detail_dictionary = {}
                    for query_species_name in query_species_list:
                        list_of_folders  = os.listdir(f"{blast_location}/{species}")
                        for folders in list_of_folders:
                            if folders.endswith(query_species_name):
                                # print()
                                
                                with open(f"{blast_location}/{species}/{folders}/{exon}/for_blast/new_query.txt", 'r') as query_file_open:
                                    # alignment_detail = query_file_open.readlines()[0]
                                    alignment_detail = query_file_open.readlines()[0].split("set")[1][:-1]
                                    alignment_detail_dictionary[query_species_name] = alignment_detail
                                with open (f"{blast_location}/{species}/{folders}/{exon}/for_alignment/{exon}_translated_genomic_sequence_{alignment_detail}.fa.hat2", 'r') as alignment_score_file:
                                    distance_score = float(alignment_score_file.readlines()[-1].rstrip())
                                try:
                                    with open (f"{blast_location}/{species}/{folders}/{exon}/for_blast/new_query_spliced.txt", 'r') as new_query_file:
                                        new_query = (new_query_file.readlines()[-1].rstrip())
                                except:
                                    new_query = "Issue with query"
                                print(query_species_name,exon, alignment_detail, distance_score,f"\n{new_query}")
                    while True:            
                        check_alignment = input("Check Alignment??")
                        if check_alignment[0].lower() == "y" or check_alignment[0].lower() == "n":
                            break
                    if check_alignment[0].lower() == "y":
                        for query_species_name in query_species_list:
                            list_of_folders  = os.listdir(f"{blast_location}/{species}")
                            for folders in list_of_folders:
                                if folders.endswith(query_species_name):
                                    with open (f"{blast_location}/{species}/{folders}/{exon}/for_alignment/alignment_clustal_{exon}_translated_genomic_sequence_{alignment_detail_dictionary[query_species_name]}.fa.txt", 'r') as clustal_alignment_file:
                                        clustal_alignment = clustal_alignment_file.readlines()
                                    print("".join(clustal_alignment))
                                    
                                # assert False
                    while True:
                        try:
                            species_number = (input(f"Choose Species 1 - {len(query_species_list)}"))
                            
                            
                            if species_number[0].lower() == 'n':
                                output +=  f'{species},{scaffold},000,000,0,Y,Error_{exon},00,00,00\n'
                                break 
                            output += coordinate_detail_dictionary[query_species_list[int(species_number)-1]][exon]
                            break
                        except:
                            print("Retry")
                    break
                elif choice[0].lower() == "n":
                    output +=  f'{species},{scaffold},000,000,0,Y,Error_{exon},00,00,00\n'
                    break
                
                
        else:
            print(exon, "Match")
            print(coordinates_list)
            output += coordinate_detail_dictionary[query_species_list[0]][exon]
        
       

        # assert False
    # print(coordinate_detail_dictionary)
    print(output)
    with open(f"{blast_location}/{species}/final_coordinates.csv", "w") as out_file:
        out_file.write(output)
    
