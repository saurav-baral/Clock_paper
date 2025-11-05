# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 09:52:58 2024

@author: sauba
"""

import os
from Bio import SeqIO
import subprocess

file_location = os.getcwd()
family_group_list = ["3.Satyrine","4.Pierinae","5.Coliadinae","6.Heliconiinae_Danainae_Nymphalinae","7.Papilionidae","8.Lycaenidae"]
for family_group in family_group_list:



    output_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/15.Timeout_Exon_Analysis/{family_group}"

    genome_location = "/mnt/griffin/saubar/Genomes_2023-12-26"



    def get_genome_file(genome_location,species):
        list_of_files_in_genome_folder = os.listdir(f"{genome_location}/{species}")
        for file in list_of_files_in_genome_folder:
            if file.endswith("_genomic.fna"):
                genome_file = file
        return(genome_file)  

    def get_gene_sequence(genome_location, 
                          species, 
                          genome_file, 
                          scaffold,
                          gene_start,
                          gene_end,
                          complement,
                          output_location):

        
        print("Getting Gene")
        list_of_folders = os.listdir(output_location)
        # print(list_of_files_inside_annotated_species_folder)
        # print(
        if f"0.Temp" not in list_of_folders:
            os.mkdir(f"{output_location}/0.Temp")

    #     subprocess.run(f'samtools faidx "{genome_location}/{species}/{genome_file}"', shell = True, stderr = subprocess.DEVNULL)
    #     print(f'samtools faidx "{genome_location}/{species}/{genome_file}" {scaffold}:{gene_start}-{gene_end}')
        subprocess.run(f'samtools faidx "{genome_location}/{species}/{genome_file}" {scaffold}:{gene_start}-{gene_end} > "{output_location}/0.Temp/temp_genome.fa"', shell = True, stderr = subprocess.DEVNULL)
        
        genome = SeqIO.parse(f"{output_location}/0.Temp/temp_genome.fa", "fasta")
        for entries in genome:
            gene_sequence = entries.seq
            if complement == "1":
                gene_sequence = gene_sequence.reverse_complement()
            break
        # print(gene_sequence)
        return (gene_sequence)
    def make_output_folder(output_location):
        list_of_folders = os.listdir(output_location)
        if "2.Final_output" in list_of_folders:
            subprocess.run(f'rm -r "{output_location}/2.Final_output"', shell = True, stderr = subprocess.DEVNULL)
        os.mkdir(f"{output_location}/2.Final_output")

    def write_exon_output(output_location,species,exon,exon_output):
        
        
        list_of_folders_in_final_output = os.listdir(f"{output_location}/2.Final_output")
        if species not in list_of_folders_in_final_output:
            os.mkdir(f"{output_location}/2.Final_output/{species}")
        
        with open(f"{output_location}/2.Final_output/{species}/{exon}.fa", 'w') as out_file:
            out_file.write(exon_output)
            
            

    make_output_folder(output_location)

    list_of_species = os.listdir(f"{output_location}/1.Blast_result")
    if 'desktop.ini' in list_of_species:
        list_of_species.remove("desktop.ini")
    # print(list_of_species)
    # list_of_species = ["Heliconius_charithonia"]
    # list_of_species = ["Maniola_hyperantus"]
    error_exons = ''
    for species in list_of_species:
        print(species)
        genome_file = get_genome_file(genome_location,species)
        print(genome_file)
        previous_coordinate_start = ''
        previous_scaffold = ''
        with open(f"{output_location}/1.Blast_result/{species}/final_coordinates.csv", 'r') as coordinate_final_file:
            coordinate_file_list = coordinate_final_file.readlines()
        last_exon = len(coordinate_file_list)-1
        
        for i in range(1,len(coordinate_file_list)):
            
            coordinate_line_split = coordinate_file_list[i].split(",")
            
                
    #         print(coodinate_line_split[5])
            if coordinate_line_split[5] == "N":
                if previous_coordinate_start == '':
                    previous_coordinate_start = coordinate_line_split[2]
                    previous_scaffold = coordinate_line_split[1]
                complement = coordinate_line_split[4]
                start,stop = coordinate_line_split[2],coordinate_line_split[3]
                print(start,stop)
                query_name = coordinate_line_split[6]
                scaffold = coordinate_line_split[1]
                exon = f"Exon_{query_name.split('_')[-1]}"
                left_oh, right_oh = query_name.split('_')[-6], query_name.split('_')[-4]
                print(left_oh, right_oh)
    #             print(complement)
                if complement == "0":
                    if start < previous_coordinate_start or scaffold != previous_scaffold:
                        print(f"Error in Coordinates\nPrevious_coordinate = {coordinate_file_list[i-1]}\nCurrent_coordinate = {coordinate_file_list[i]}")
                        prompt = input("Continue (Y/N)?")
                        if prompt[0].lower() == "n":
                            assert False
                elif complement == "1":
                    if start > previous_coordinate_start or scaffold != previous_scaffold:
                        print(f"Error in Coordinates\nPrevious_coordinate = {coordinate_file_list[i-1]}\nCurrent_coordinate = {coordinate_file_list[i]}")
                        prompt = input("Continue (Y/N)?")
                        if prompt[0].lower() == "n":
                            assert False
                else:
                    print("ERROR")
                
            elif coordinate_line_split[5] == "Y":
                query_name = coordinate_line_split[6]
                exon = f"Exon_{query_name.split('_')[-1]}"
                scaffold = coordinate_line_split[1]
                error_exons += f"{species},{scaffold},000,000,0,Y,Error_{exon},00,00,00\n"
                continue
            print(coordinate_file_list[i])  
            gene_sequence = get_gene_sequence(genome_location, 
                                              species, 
                                              genome_file, 
                                              scaffold,
                                              start,
                                              stop,
                                              complement,
                                              output_location)
            
            new_start = int(start) - 2
            new_stop = int(stop) +2
            # elif complement == "1":
            #     new_start = int(start) +2
            #     new_stop = int(stop) -2
            # else:
            #     print("ERROR")
            #     assert False
            
            gene_sequence_to_check = get_gene_sequence(genome_location, 
                                              species, 
                                              genome_file, 
                                              scaffold,
                                              new_start,
                                              new_stop,
                                              complement,
                                              output_location)
            # print(gene_sequence, "\n\n\n", gene_sequence_to_check, "\n\n\n")
            
            # print(gene_sequence_to_check [:2].upper(), gene_sequence_to_check [-2:].upper())
            # assert False
            current_exon = query_name.split('_')[-1]
            print(f"Current Exon = {current_exon}")
            if current_exon == "1":
                if complement == 1:
                    translated = gene_sequence.reverse_complement.translate()
                else:
                    translated = gene_sequence.translate()
                print(translated)
                if "*" in translated and translated[0] != "M" and gene_sequence_to_check [-2:].upper() != "GT":
                    error_exons += f"{species},{scaffold},000,000,0,Y,Error_{exon},00,00,00\n"
                else:
                    exon_output = f">{species}_{exon}_{scaffold}_{start}_{stop}_left_{left_oh}_right_{right_oh}\n{gene_sequence}"
                    print(exon_output)
            
                    write_exon_output(output_location,species,exon,exon_output)
            elif current_exon == f"{last_exon}":
                if complement == 1:
                    translated = gene_sequence[:-int(left_oh)].reverse_complement.translate()
                    print(gene_sequence)
                else:
                    translated = gene_sequence[int(left_oh):].translate()
                    print(gene_sequence)

                if "*" in translated[:-1] and translated[-1] != "*" and gene_sequence_to_check [:2].upper() != "AG":
                    error_exons += f"{species},{scaffold},000,000,0,Y,Error_{exon},00,00,00\n"
                    print(translated)
                else:
                    
                    exon_output = f">{species}_{exon}_{scaffold}_{start}_{stop}_left_{left_oh}_right_{right_oh}\n{gene_sequence}"
                    print(exon_output)
                    write_exon_output(output_location,species,exon,exon_output)
    #         assert False
            else:
                print(gene_sequence_to_check [:2].upper(), gene_sequence_to_check [-2:].upper())
                if gene_sequence_to_check [:2].upper() != "AG" and gene_sequence_to_check [-2:].upper() != "GT":
                    
                    error_exons += f"{species},{scaffold},000,000,0,Y,Error_{exon},00,00,00\n"
                else:
                    exon_output = f">{species}_{exon}_{scaffold}_{start}_{stop}_left_{left_oh}_right_{right_oh}\n{gene_sequence}"
                    print(exon_output)
        
                    write_exon_output(output_location,species,exon,exon_output)
                
    #         assert False
            previous_coordinate_start = start
            previous_scaffold = scaffold
        error_exons += "\n\n"
    #     print(coordinate_file_list)
    #     break
    coordinate_file = ''
    print(error_exons)
    with open(f"{output_location}/error_exons.txt", 'w') as error_out_file:
        error_out_file.write(error_exons)