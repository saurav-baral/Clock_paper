# -*- coding: utf-8 -*-
"""
Created on Wed May 29 13:13:58 2024

@author: sauba
"""






import os
from Bio import SeqIO
import io
from Bio.Seq import Seq
import subprocess
import sys


file_location = os.getcwd()
family_group = file_location.split("/")[-1]


location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/13.Cryptochrome_1b_Exon_Analysis/{family_group}"
genome_location = "/mnt/griffin/saubar/Genomes_2023-12-26"


species_list = [f"{sys.argv[1]}"]
blast_output_location = f"{location}/1.Blast_result"
annotated_genome_location = blast_output_location
query_location = f"{location}/1.Query"



def make_raw_files_for_alignment(gene_sequence,annotated_genome_location,annotated_species_name,error_exon,query_fasta_sequence,query_length):
    for offset in range(3):
        translated_sequence = str(gene_sequence[offset:].translate()).split("*")
        for i in range(len(translated_sequence)):
            if len(translated_sequence[i])> 0.8*query_length and "XXXXXX" not in translated_sequence[i]:
                sequence_set = f">set{i+1}_frame{offset}\n{translated_sequence[i]}\n\n"
                # print(i+1, offset)
                
                with open(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_alignment/{error_exon}_translated_genomic_sequence_{i+1}_frame{offset}.fa",'w') as out_file:
                    output = f"{query_fasta_sequence}\n\n{sequence_set}"
                    out_file.write(output)
                    
def run_mafft(annotated_genome_location,annotated_species_name,error_exon):
    location = f'{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_alignment'
    list_of_files_to_run_mafft_on = os.listdir(location)
    for file in list_of_files_to_run_mafft_on:
        if file.endswith(".fa"):
            command = f'"mafft" --localpair --maxiterate 16 --reorder --distout "{location}/{file}" > "{location}/alignment_{file}.txt"'
            # print(command)
            subprocess.run(f'{command}', shell=True, stderr = subprocess.DEVNULL) 
            # os.system(f'{command}')
            command = f'"mafft" --localpair --clustalout --maxiterate 16 --reorder "{location}/{file}" > "{location}/alignment_clustal_{file}.txt"'
            subprocess.run(f'{command}', shell=True, stderr = subprocess.DEVNULL) 
            
            # os.system(f'{command}')
    return(location)

def process_mafft_output(mafft_run_folder, error_exon):
    list_of_files_in_mafft_run_folder = os.listdir(mafft_run_folder)
    score_output = []
    score = 99
    min_score_sequence = ''
    alignment_file = ''
    for file in list_of_files_in_mafft_run_folder:
        if file.endswith(".fa.hat2"):
            with io.open(f"{mafft_run_folder}/{file}", 'r') as dist_matrix_file:
                dist_matrix_list = dist_matrix_file.readlines()
            
            sequence_name = dist_matrix_list[-2].rstrip().split("=")[1]
            distance_score = float(dist_matrix_list[-1].rstrip())
            if len(score_output) < 5:
                score_output.append( [sequence_name,distance_score])
            else:
                for i in range(len(score_output)):
                    score_at_this_index = score_output[i][1]
                    if distance_score < score_at_this_index:
                        score_output[i] = [sequence_name,distance_score]
                        break
            if distance_score < score:
                score = distance_score
                min_score_sequence = sequence_name
                alignment_file = f'alignment_{file.replace(".hat2",".txt")}'
    # print(score_output)
    print(f"min = {min_score_sequence}, {score}" )
    # score_out_merged = '\n'.join(score_output)
    print(f"5 top scores:\n{score_output}")
    print(f"{mafft_run_folder}/{alignment_file}")
    
    alignment_file = SeqIO.parse(f"{mafft_run_folder}/{alignment_file}", 'fasta')
    # print (records.id)
    fasta_start_position = 0
    fasta_end_position = 0
    start_switch = 0
    end_switch = 0
    alignment_name = ''
    for records in alignment_file:
        
        
        
        if error_exon in records.id:
            gap_counter = 0
            base_counter = 0
            for current_position in range(len(records.seq)):
                sequence_length = len(records.seq) - records.seq.count('-')
                # print(f"fasta_end_position {fasta_end_position} fasta_start_position {fasta_start_position}")
                # print(f"current_position = {current_position}, {len(records.seq)}")
                # print(records.seq[current_position])
                # print("fasta_start_position",fasta_start_position)
                # print("start_switch",start_switch)
                
                # print(gap_counter, base_counter)
                # print(5,0.2*sequence_length)
                if start_switch == 1 and records.seq[current_position] == "-" and gap_counter > 3 and (len(records.seq[:current_position]) - records.seq[:current_position].count('-')) < (0.1*sequence_length):
                #  :
                    # print("\n\nhere\n\n")
                    start_switch = 0
                    gap_counter = 0

                if "-" not in records.seq[current_position] and start_switch == 0 :
                    fasta_start_position = current_position
                    start_switch = 1
                
                    
                
                    
                if end_switch == 1 and "-" not in records.seq[current_position] and (len(records.seq[current_position:]) - records.seq[current_position:].count('-')) > 0.1*sequence_length :
                    end_switch = 0
                    gap_counter = 0
                
                    
                if start_switch == 1 and records.seq[current_position] == "-" and end_switch == 0:
                    # print(f"base_counter {base_counter}")
                    # print(f"fasta_end_position {fasta_end_position}")
                    
                    # print(f"fasta_end_position {fasta_end_position}")
                    
                    fasta_end_position = current_position
                    end_switch = 1
                if "-" in records.seq[current_position]:
                    gap_counter += 1
                else:
                    gap_counter = 0
                    base_counter += 1
        else:
            alignment_name = records.id
            if fasta_end_position == 0:
                
                fasta_end_position = current_position
            end_switch = 1
            start_switch = 1

        if start_switch == 1 and end_switch == 1:
            print(fasta_start_position, fasta_end_position)
            gene_sequence = records.seq[fasta_start_position:fasta_end_position]
            print(f"{records.id}\n{gene_sequence}")
    return(gene_sequence, alignment_name)

def get_genome_file(genome_location,species):
    list_of_files_in_genome_folder = os.listdir(f"{genome_location}/{species}")
    for file in list_of_files_in_genome_folder:
        if file.endswith("_genomic.fna"):
            genome_file = file
    return(genome_file)  

def get_annotated_genome_name(annotated_genome_location, species):
    list_of_annotated_genomes = os.listdir(annotated_genome_location)
  
    # print(annotated_species)
    for annotated_species in list_of_annotated_genomes:
        if annotated_species.endswith(species):
            return(annotated_species)
             
    if annotated_species_name == '':
        print(f"Error with annotated species name")
        assert False

def check_and_make_folders(annotated_genome_location,annotated_species_name,error_exon):
    list_of_files_inside_annotated_species_folder = os.listdir(f"{annotated_genome_location}/{annotated_species_name}")

    if f"Period_gene_genomic_sequence_individual_exon_{query_species}" not in list_of_files_inside_annotated_species_folder:
        os.mkdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}")
    
    list_of_exon_directories = os.listdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/")
    
    if error_exon not in list_of_exon_directories:
        os.mkdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}")
        os.mkdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_alignment")
    elif error_exon in list_of_exon_directories:
        list_of_folders_1 = os.listdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}")
        if "for_alignment" not in list_of_folders_1:
            os.mkdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_alignment")
        list_of_files = os.listdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_alignment")
        for file in list_of_files:
            os.remove(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_alignment/{file}")

def get_gene_sequence(genome_location, species, genome_file, scaffold,gene_start,gene_end,complement,annotated_genome_location,annotated_species_name):
    from Bio import SeqIO
    print("Getting Gene")
    list_of_files_inside_annotated_species_folder = os.listdir(f"{annotated_genome_location}/{annotated_species_name}")
    # print(list_of_files_inside_annotated_species_folder)
    # print(
    if f"Period_gene_genomic_sequence_individual_exon_{query_species}" not in list_of_files_inside_annotated_species_folder:
        os.mkdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}")

    list_of_files_inside_indiv_exon_folder = os.listdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}")
    if "temp" not in list_of_files_inside_indiv_exon_folder:
         os.mkdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/temp")

    # os.system(f'samtools faidx "{genome_location}/{species}/{genome_file}"')
    # subprocess.run("pwd")
    subprocess.run(f'samtools faidx "{genome_location}/{species}/{genome_file}"', shell = True, stderr = subprocess.DEVNULL)
    # os.system(f'samtools faidx "{genome_location}/{species}/{genome_file}" {scaffold}:{gene_start}-{gene_end} > "{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/temp/temp_genome.fa"')
    subprocess.run(f'samtools faidx "{genome_location}/{species}/{genome_file}" {scaffold}:{gene_start}-{gene_end} > "{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/temp/temp_genome.fa"', shell = True, stderr = subprocess.DEVNULL)
    genome = SeqIO.parse(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/temp/temp_genome.fa", "fasta")
    for entries in genome:
        gene_sequence = entries.seq
        if complement == "1":
            gene_sequence = gene_sequence.reverse_complement()
        break
    # print(gene_sequence)
    return (gene_sequence)


def mafft_process(previous_exon_coordinates,
                  next_exon_coordinates,
                  current_exon_coordinates,
                  query_species,
                  query_transcript,
                  query_location,
                  annotated_genome_location,
                 genome_location):
    query_exon = current_exon_coordinates[6].split("query")[-1][1:]
    error_exon = query_exon
    # print("AA")
    genome_file = get_genome_file(genome_location,species)
    
    upstream_exon_line, downstream_exon_line = previous_exon_coordinates.split(","),next_exon_coordinates.split(",")
    complement,scaffold = upstream_exon_line[4], upstream_exon_line[1]

    if complement == "0":
        gene_start = min(int(upstream_exon_line[2]),int(upstream_exon_line[3]))
        gene_end = max(int(downstream_exon_line[2]),int(downstream_exon_line[3]))
    
    if complement == "1":
        gene_start = min(int(downstream_exon_line[2]),int(downstream_exon_line[3]))
        gene_end = max(int(upstream_exon_line[2]),int(upstream_exon_line[3]))
    
    
    with open(f"{query_location}/{query_species}/{query_transcript}/query_{query_exon}.fa", 'r') as query_file:
        query_file_list = query_file.readlines()
        query_fasta_sequence = "".join(query_file_list)
        left_overhang = query_file_list[0].split("Frame")[1][1]
        right_overhang = query_file_list[0].split("rightoh")[1][1]
        original_query_name = query_file_list[0]
        original_query = query_file_list[1]
        
    query_length = len(query_fasta_sequence.split("\n")[1])
    
    annotated_species_name = get_annotated_genome_name(annotated_genome_location, species)
    # print(annotated_species_name)

    # print(annotated_genome_location)
    gene_sequence = get_gene_sequence(genome_location, species, genome_file, scaffold,gene_start,gene_end,complement,annotated_genome_location,annotated_species_name)
    # assert False
    check_and_make_folders(annotated_genome_location,annotated_species_name,error_exon)
    
    
    make_raw_files_for_alignment(gene_sequence,annotated_genome_location,annotated_species_name,error_exon,query_fasta_sequence,query_length)
    
    return(annotated_species_name, error_exon,species, genome_file,gene_sequence,left_overhang,right_overhang,scaffold, original_query_name,original_query )




if "desktop.ini" in species_list:
    species_list.remove("desktop.ini")

for species in species_list:
    list_of_query_species = os.listdir(query_location)
    if "desktop.ini" in list_of_query_species:
        list_of_query_species.remove("desktop.ini")
#     list_of_query_species = ["Papilio_xuthus"]
    for query_species in list_of_query_species:
        list_of_query_transcripts = os.listdir(f"{query_location}/{query_species}")
        if "desktop.ini" in list_of_query_transcripts:
            list_of_query_transcripts.remove("desktop.ini")
        for query_transcript in list_of_query_transcripts:


            



            output_coordinate_file = "Species," + "Scaffold," + "Start," + "Stop," + "Complement," + "Error," + "Gene,"+ "Query_start," + "Query_stop,"+ "Query_Length," +  "AG_GT," + "Spliceator_prediction\n"
            #exons_to_check_list = ["Exon_2","Exon_3","Exon_4","Exon_5","Exon_6","Exon_7","Exon_8","Exon_9","Exon_10","Exon_11","Exon_12"]
#             exons_to_check_list = ["Exon_9"]
            exons_to_check_list = []
            list_of_files_in_species_folder = os.listdir(f"{blast_output_location}/{species}")
            coordinate_file_name = ''
        #     print(list_of_files_in_species_folder)
            for file_names in list_of_files_in_species_folder:

                if file_names.endswith("_coordinates_old.csv"):
                    coordinate_file_name = file_names
            if coordinate_file_name =='':
                print(f"Coordinate file error")
                assert False

            with io.open(f"{blast_output_location}/{species}/{coordinate_file_name}", 'r') as temp_file_open:
                coordinate_file_lines = temp_file_open.readlines()
            # print(coordinate_file_lines[0])
            for i in range(2,len(coordinate_file_lines)-1):
                exons_to_check_list.append(f"Exon_{i}")
            print(exons_to_check_list)
            
            for exons_to_check in exons_to_check_list:

                for i in range(1,len(coordinate_file_lines)):
                    coordinate_file_lines_split = coordinate_file_lines[i].split(",")
                    current_exon_coordinates = coordinate_file_lines[i].split(",")
                    if i == 1:
                        scaffold_first_exon = coordinate_file_lines[i].split(",")[1]
                        first_exon_start, first_exon_end = coordinate_file_lines[i].split(",")[2],coordinate_file_lines[i].split(",")[3]
                    # if coordinate_file_lines_split[5] == "Y":
                    if  coordinate_file_lines_split[6].endswith(exons_to_check):
                        if i == 1:
                            print(f'First Exon has errors\n{coordinate_file_lines[i]}')
                            continue
#                             assert False
                        else:
                            print(f"Scaffold = {scaffold_first_exon}, start = {first_exon_start}, end = {first_exon_end}")
                            print(coordinate_file_lines[i])
                            # process_current_exon = input("Error Found! Process?")
                            process_current_exon = 'y'
                            if process_current_exon.lower()[0] == "y":
                                if i>2:
                                    previous_exon_number = i-2
                                else:
                                    previous_exon_number = i-1
                                if i+2 <=int(exons_to_check_list[-1].split("_")[1]):
                                    next_exon_number = i+2
                                else:
                                    next_exon_number = i+1
                                if i >= int(exons_to_check_list[-1].split("_")[1]):
                                    next_exon_number = int(exons_to_check_list[-1].split("_")[1])
#                                     continue

                                '''
                                the above code ensures that missing exons do not cause issues
                                '''

                                previous_exon_coordinates = ''
                                next_exon_coordinates = ''
                                while True:
                                    print(previous_exon_number)
                                    if coordinate_file_lines[previous_exon_number].split(",")[5] == "N":
                                        previous_exon_coordinates = coordinate_file_lines[previous_exon_number]
                                        break
                                    # else:
                                    #     proceed_prompt = input(f"Previous exon coordinate\n{coordinate_file_lines[previous_exon_number]}\nProceed?")
                                    #     if proceed_prompt.lower()[0] == "y":
                                    #         previous_exon_coordinates = coordinate_file_lines[previous_exon_number]
                                    #         break
                                    previous_exon_number -= 1
                                while True:
                                    print(next_exon_number)
                                    try:
                                        if coordinate_file_lines[next_exon_number].split(",")[5] == "N":
                                            next_exon_coordinates = coordinate_file_lines[next_exon_number]
                                            break
                                    except:
                                        if coordinate_file_lines[next_exon_number-1].split(",")[5] == "N":
                                            next_exon_coordinates = coordinate_file_lines[next_exon_number-1]
                                            break
                                   
                                    next_exon_number += 1


                                print(previous_exon_coordinates)
                                print(next_exon_coordinates)
                                annotated_species_name, error_exon, species, genome_file, gene_sequence,left_overhang,right_overhang,scaffold,original_query_name, original_query = mafft_process(previous_exon_coordinates,
                                                                                                                                  next_exon_coordinates,
                                                                                                                                  current_exon_coordinates,
                                                                                                                                  query_species,
                                                                                                                                  query_transcript,
                                                                                                                                  query_location,
                                                                                                                                  annotated_genome_location,
                                                                                                                                     genome_location)

                                mafft_run_folder = run_mafft(annotated_genome_location,annotated_species_name,error_exon)

                                possible_gene_sequence, alignment_name = process_mafft_output(mafft_run_folder, error_exon)

                                print(f"Original query = {original_query}")
                                print(f"New query = {possible_gene_sequence}")
                                
                                list_of_folders_inside_exon_folder = os.listdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}")
                                if "for_blast" not in list_of_folders_inside_exon_folder:
                                    os.mkdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_blast")
                                else:
                                    list_of_files = os.listdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_blast")
                                    for file in list_of_files:
                                        os.remove(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_blast/{file}")


                                with open(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_blast/new_query.txt",'w') as query_file:
                                    query = f">Query_{error_exon}_{alignment_name}\n{possible_gene_sequence}"
                                    query_file.write(query)
