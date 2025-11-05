# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 12:47:44 2024

@author: sauba
"""
import os
import subprocess

file_location = os.getcwd()
family_group = file_location.split("/")[-1]

location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/14.Cryptochrome_2_Exon_Analysis/{family_group}"

species_list = os.listdir(f"{location}/1.Blast_result")
if "desktop.ini" in species_list:
    species_list.remove("desktop.ini")
error_species = ''
# species_list=["Pararge_aegeria"]



for species in species_list:
    
    list_of_files_in_species_result_folder = os.listdir(f"{location}/1.Blast_result/{species}")
    backup_present = 0
    for files in list_of_files_in_species_result_folder:
        if files.endswith("_coordinates_old_backup.csv"):
            backup_present = 1

    if backup_present == 0:
        subprocess.run(f'cp "{location}/1.Blast_result/{species}/{species}_coordinates_old.csv" "{location}/1.Blast_result/{species}/{species}_coordinates_old_backup.csv"', shell = True)

    
    with open(f"{location}/1.Blast_result/{species}/{species}_coordinates_old_backup.csv", 'r') as coordinate_file:
        coordinate_list = coordinate_file.readlines()
    print("".join(coordinate_list))
    input("continue?")
    first_exon_error = coordinate_list[1].split(",")[5]
    last_exon_error = coordinate_list[-1].split(",")[5]
    
    
    
        
    if first_exon_error == "Y":
         
        print(f"\n\nFirst Exon: \n{coordinate_list[1]}")
        exon_name = coordinate_list[1].split(",")[6]
        while (True):

            fix_option = input("Fix!!\n1.Add New Coordinate\n2.Just change error from 'Y' to 'N'\n3.Stop 'n'")

            if fix_option == "n":
                assert False
            elif fix_option == "1":

                for i in range(len(coordinate_list)):
                    if coordinate_list[i].split(",")[5] == "N":
                        Species,Scaffold,Start,Stop,Complement = coordinate_list[i].split(",")[:5]
                        print(Species,Scaffold,Start,Stop,Complement)
                        if Complement == "0":
                            new_start = int(Start) - 2000*i
                            new_stop = int(Stop) - 2000*i
                        elif Complement == "1":
                            new_start = int(Start) + 2000*i
                            new_stop = int(Stop) + 2000*i

                        coordinate_list[1] = f"{Species},{Scaffold},{new_start},{new_stop},{Complement},N,{exon_name},0,0,0\n"             
                        break
                error_species +=f"{Species},{Scaffold},000,000,0,Y,Error_Exon_1,00,00,00\n"    
                break
            

            elif fix_option == "2":
                print("changing Y to N")
                
                coordinate_list[1] = coordinate_list[1].replace(',Y,',',N,')
                
                break
            else:
                print("Retry?")
        

    if last_exon_error == "Y":
        print(f"\n\nLast Exon: \n{coordinate_list[-1]}")
        exon_name = coordinate_list[-1].split(",")[6]
        exon_number = exon_name.split("_")[-1]
        while (True):

            fix_option = input("Fix!!\n1.Add New Coordinate\n2.Just change error from 'Y' to 'N'\n3.Stop 'n'")

            if fix_option == "n":
                assert False
            elif fix_option == "1":

                for i in range(len(coordinate_list)):
                    position = len(coordinate_list)-1-i
                    if coordinate_list[position].split(",")[5] == "N":
                        Species,Scaffold,Start,Stop,Complement = coordinate_list[position].split(",")[:5]
    #                     print(Species,Scaffold,Start,Stop,Complement)
                        if Complement == "0":
                            new_start = int(Start) + 2000*i
                            new_stop = int(Stop) + 2000*i
                        elif Complement == "1":
                            new_start = int(Start) - 2000*i
                            new_stop = int(Stop) - 2000*i

                        coordinate_list[-1] = f"{Species},{Scaffold},{new_start},{new_stop},{Complement},N,{exon_name},0,0,0\n"             
                        
                        break
                error_species +=f"{Species},{Scaffold},000,000,0,Y,Error_Exon_{exon_number},00,00,00\n"
                break

            elif fix_option == "2":
                coordinate_list[-1] = coordinate_list[-1].replace(',Y,',',N,')
                break
            else:
                print("Retry?")
        

    new_coordinates = "".join(coordinate_list)
    print(new_coordinates)

    with open(f"{location}/1.Blast_result/{species}/{species}_coordinates_old.csv", 'w') as out_coordinate_file:
        out_coordinate_file.write(new_coordinates)
    
with open(f"{location}/error_in_first_and_last.txt", 'w') as error_file:
    error_file.write(error_species)
print(error_species)