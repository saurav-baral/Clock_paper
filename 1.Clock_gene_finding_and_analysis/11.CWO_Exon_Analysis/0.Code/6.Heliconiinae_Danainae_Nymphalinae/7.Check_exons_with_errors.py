# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 12:41:19 2024

@author: sauba
"""

import os
from Bio import SeqIO
import subprocess

file_location = os.getcwd()
family_group = file_location.split("/")[-1]



output_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/11.CWO_Exon_Analysis/{family_group}"
genome_location = "/mnt/griffin/saubar/Genomes_2023-12-26"

def make_folder_error_exon(output_location,species, error_exon):
    import subprocess
    list_of_folders = os.listdir(f"{output_location}/1.Blast_result/{species}")
    if "Error_exon_processing" not in list_of_folders:
        os.mkdir(f"{output_location}/1.Blast_result/{species}/Error_exon_processing")
        
    list_of_folder_2 = os.listdir(f"{output_location}/1.Blast_result/{species}/Error_exon_processing")
    if error_exon in list_of_folder_2:
        subprocess.run(f'rm -r "{output_location}/1.Blast_result/{species}/Error_exon_processing/{error_exon}"', shell = True, stderr = subprocess.DEVNULL)
    os.mkdir(f"{output_location}/1.Blast_result/{species}/Error_exon_processing/{error_exon}")
    
    list_of_folders_3 = os.listdir(f"{output_location}/1.Blast_result/{species}/Error_exon_processing/{error_exon}")
    if "Temp_query" not in list_of_folders_3:
        os.mkdir(f"{output_location}/1.Blast_result/{species}/Error_exon_processing/{error_exon}/Temp_query")
        
def make_query(output_location,error_exon,species):
    from Bio import SeqIO
    import random
    import os
    exon = error_exon.split("Error")[1][1:]
    print(exon)
#     print("Here Here")
    initial_list_of_species = os.listdir(f"{output_location}/2.Final_output")
    
    
    try:
        initial_list_of_species.remove("desktop.ini")
    except:
        pass
    list_of_species = []
#     print(initial_list_of_species)
    for species_name in initial_list_of_species:
        list_of_files_in_species_folder = os.listdir(f"{output_location}/2.Final_output/{species_name}")
#         print(list_of_files_in_species_folder)
#         print(f"{exon}.fa")
        if f"{exon}.fa" in os.listdir(f"{output_location}/2.Final_output/{species_name}"):
            list_of_species.append(species_name)
    
    
    if len(list_of_species) > 6:
        list_of_species = random.sample(list_of_species, 5)
    
#    add_species = input("Preferred Species?")
#    if add_species.lower() != "n":
#        if add_species not in list_of_species:
#            list_of_species.append(add_species)
#     list_of_species.append("Helleia_helle")
    for query_species in list_of_species:
        print("current_query :", query_species)
        list_of_exons = os.listdir(f"{output_location}/2.Final_output/{query_species}")
        try:
            list_of_exons.remove("desktop.ini")
        except:
            pass
#         print(query_species,list_of_exons)
        if f"{exon}.fa" in list_of_exons:
            print (f"{exon}.fa")
            exon_file = SeqIO.parse(f"{output_location}/2.Final_output/{query_species}/{exon}.fa", 'fasta')
            for records in exon_file:
                print(records.id)
                print(records.seq[int(records.id.split("_")[-3]):].translate())
                left_oh = records.id.split("_")[-3]
                right_oh = records.id.split("_")[-1]
    #             print(left_oh)
                query_sequence = records.seq[int(records.id.split("_")[-3]):].translate()
            with open(f"{output_location}/1.Blast_result/{species}/Error_exon_processing/{error_exon}/Temp_query/{query_species}_{exon}.fa", "w") as out_file:
                output = f">{records.id}\n{query_sequence}"
                out_file.write(output)
    return(left_oh,right_oh)

def get_genomic_coordinates(output_location,error_exon,species):
    error_exon_number = error_exon.split("_")[-1]
#     print(error_exon_number)
#     return
    with open(f"{output_location}/1.Blast_result/{species}/final_coordinates.csv", 'r') as final_coordinate_file:
        final_coordinate_file_lines = final_coordinate_file.readlines()
    for i in range(1, len(final_coordinate_file_lines)):
#         print(final_coordinate_file_lines[i])
        current_position = i
        next_position = i
#         print(final_coordinate_file_lines[i].split(",")[6].split("Exon")[1][1:].strip())
#         print(error_exon_number)
        scaffold = final_coordinate_file_lines[current_position-1].split(",")[1]

        original_query_name =f"{'_'.join(final_coordinate_file_lines[current_position-1].split(',')[6].split('_')[:-1])}_{current_position}"
        
        
        if final_coordinate_file_lines[i].split(",")[6].split("Exon")[1][1:].strip() == error_exon_number:
#             print(current_position, next_position)
    
            if error_exon_number == "1":
                current_position = current_position
                while final_coordinate_file_lines[next_position+1].split(",")[5] == "Y":
                        next_position += 1
                complement = final_coordinate_file_lines[current_position].split(",")[4]
                
                if complement == "0":
                    fragment_start = int(final_coordinate_file_lines[current_position].split(",")[3])-1000
                    fragment_end = final_coordinate_file_lines[next_position+1].split(",")[2]
                if complement == "1":
                    fragment_start = final_coordinate_file_lines[next_position+1].split(",")[3]
                    fragment_end = int(final_coordinate_file_lines[current_position].split(",")[2]) + 1000

                scaffold = final_coordinate_file_lines[current_position].split(",")[1]
                
                
                
            elif error_exon_number == "16":
                next_position = next_position
                while final_coordinate_file_lines[current_position-1].split(",")[5] == "Y":
                    current_position = current_position-1
                complement = final_coordinate_file_lines[current_position].split(",")[4]
                
                if complement == "0":
                    fragment_start = int(final_coordinate_file_lines[current_position-1].split(",")[3])
                    fragment_end = int(final_coordinate_file_lines[next_position].split(",")[2])+1000
                if complement == "1":
                    fragment_start = int(final_coordinate_file_lines[next_position].split(",")[3])-1000
                    fragment_end = final_coordinate_file_lines[current_position-1].split(",")[2] 
                
                scaffold = final_coordinate_file_lines[current_position].split(",")[1]
                original_query_name =f"{'_'.join(final_coordinate_file_lines[current_position].split(',')[6].split('_')[:-1])}_{current_position}"
                
            else: 
                while final_coordinate_file_lines[current_position-1].split(",")[5] == "Y":
                    current_position = current_position-1
                while final_coordinate_file_lines[next_position+1].split(",")[5] == "Y":
                    next_position += 1

                complement = final_coordinate_file_lines[current_position-1].split(",")[4]
                original_query_name =f"{'_'.join(final_coordinate_file_lines[current_position-1].split(',')[6].split('_')[:-1])}_{current_position}"

                if complement == "0":
                    fragment_start = final_coordinate_file_lines[current_position-1].split(",")[3]
                    fragment_end = final_coordinate_file_lines[next_position+1].split(",")[2]
                if complement == "1":
                    fragment_start = final_coordinate_file_lines[next_position+1].split(",")[3]
                    fragment_end = final_coordinate_file_lines[current_position-1].split(",")[2]
            
                scaffold = final_coordinate_file_lines[current_position-1].split(",")[1]

                original_query_name =f"{'_'.join(final_coordinate_file_lines[current_position-1].split(',')[6].split('_')[:-1])}_{current_position}"
                
    #             query_name =f"{''.join(final_coordinate_file_lines[current_position-1].split(''))}"
    #             print(query_name)                       
                                   
            return(fragment_start, fragment_end, scaffold,complement,original_query_name)
        

def get_genome_file(genome_location,species):
    list_of_files_in_genome_folder = os.listdir(f"{genome_location}/{species}")
    for file in list_of_files_in_genome_folder:
        if file.endswith("_genomic.fna"):
            genome_file = file
    return(genome_file) 

def get_gene_sequence(genome_location, 
                      species,                       
                      scaffold,
                      gene_start,
                      gene_end,
                      complement,
                      output_location):

    
    print("Getting Gene")
    
    genome_file = get_genome_file(genome_location,species)
    print(genome_file)
    list_of_folders = os.listdir(output_location)
    print(list_of_folders)
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

def make_raw_files_for_alignment(gene_sequence,output_location,species,error_exon):
    
    
    error_exon_location = f"{output_location}/1.Blast_result/{species}/Error_exon_processing/{error_exon}"
    list_of_folders = os.listdir(error_exon_location)
    if "for_alignment" not in list_of_folders:
        
        os.mkdir(f"{error_exon_location}/for_alignment")
    
    list_of_files = os.listdir(f"{error_exon_location}/for_alignment")
    for file in list_of_files:
        os.remove(f"{error_exon_location}/for_alignment/{file}")

    
    query_location = f"{error_exon_location}/Temp_query/"
    list_of_queries = os.listdir(query_location)
    print(list_of_queries)
#     return
    query_fasta_sequence = ''
    query_length = 0
    
    for query in list_of_queries:
        
        query_fasta_file = SeqIO.parse(f"{query_location}/{query}", 'fasta')
        for records in query_fasta_file:
            
            query_fasta_sequence = f">{records.id}\n{records.seq}\n\n"
            query_length_new = len(records.seq)
            if query_length == 0:
                query_length = len(records.seq)
            if query_length_new < query_length:
                query_length = query_length_new
                
            query_species = "_".join(records.id.split("_")[:2])
#             print(query_species)
#             return
            print(error_exon)
        #     return

        #     query_length = query_length/3
            for offset in range(3):
                translated_sequence = str(gene_sequence[offset:].translate()).split("*")
                for i in range(len(translated_sequence)):

                    if len(translated_sequence[i])> 0.8*query_length and  str(translated_sequence[i]).count("X") < 5:
                        sequence_set = f">set{i+1}_frame{offset}_{query_species}\n{translated_sequence[i]}\n\n"
        #                 print(i+1, offset)

                        with open(f"{error_exon_location}/for_alignment/{error_exon}_query_{query_species}_translated_genomic_sequence_{i+1}_frame{offset}.fa",'w') as out_file:
                            output = f"{query_fasta_sequence}\n\n{sequence_set}"
                            out_file.write(output)

def run_mafft(output_location,species,error_exon):
    error_exon_location = f"{output_location}/1.Blast_result/{species}/Error_exon_processing/{error_exon}"
    location = f'{error_exon_location}/for_alignment/'
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

def process_mafft_output(mafft_run_location, error_exon):
    error_exon = error_exon.split("Error")[1][1:]
    list_of_files_in_mafft_run_folder = os.listdir(mafft_run_location)
    score_output = []
    score = 99
    min_score_sequence = ''
    alignment_file = ''
    for file in list_of_files_in_mafft_run_folder:
        if file.endswith(".fa.hat2"):
            with open(f"{mafft_run_location}/{file}", 'r') as dist_matrix_file:
                dist_matrix_list = dist_matrix_file.readlines()
            
            sequence_name = f'set{file.split("_")[-2]}_{file.split("_")[-1].split(".")[0]}'
            query_species_name =f'{file.split("_")[4]}_{file.split("_")[5]}'
            distance_score = float(dist_matrix_list[2].strip())
            if len(score_output) < 10:
                alignment_file_now = f'alignment_{file.replace(".hat2",".txt")}'
                score_output.append( [sequence_name,distance_score,alignment_file_now])
            else:
                for i in range(len(score_output)):
                    score_at_this_index = score_output[i][1]
                    if distance_score < score_at_this_index:
                        alignment_file_now = f'alignment_{file.replace(".hat2",".txt")}'
                        score_output[i] = [sequence_name,distance_score,alignment_file_now]
                        break
            if distance_score < score:
                score = distance_score
                min_score_sequence = sequence_name
                alignment_file = f'alignment_{file.replace(".hat2",".txt")}'
#     print(score_output)
    print(f"min = {min_score_sequence}, {score}" )
    
    # score_out_merged = '\n'.join(score_output)
    print(f"5 top scores:\n{score_output}")
    
    for i in range(len(score_output)):
        temp_align_file = score_output[i][2]
        clustal_alignment_file = temp_align_file.replace("alignment_","alignment_clustal_")
        with open(f"{mafft_run_location}{clustal_alignment_file}", 'r') as clustal_file_open:
            print("".join(clustal_file_open.readlines()))
    
    print(f"{mafft_run_location}/{alignment_file}")
#     return
    
    print(f"Alignmnet file: {alignment_file}")
    clustal_alignment_file = alignment_file.replace("alignment_","alignment_clustal_")
    print(clustal_alignment_file)
    alignment_file = SeqIO.parse(f"{mafft_run_location}/{alignment_file}", 'fasta')
    # print (records.id)
    fasta_start_position = 0
    fasta_end_position = 0
    start_switch = 0
    end_switch = 0
    alignment_name = ''
    for records in alignment_file:
        
        if start_switch == 0 and end_switch == 0:
            print(error_exon)
            if error_exon in records.id:
                print(records.id)
                gap_counter = 0
                base_counter = 0
                for current_position in range(len(records.seq)):
                    sequence_length = len(records.seq) - records.seq.count('-')
#                     print(f"fasta_end_position {fasta_end_position} fasta_start_position {fasta_start_position}")
                    # print(f"current_position = {current_position}, {len(records.seq)}")
                    # print(records.seq[current_position])
                    # print("fasta_start_position",fasta_start_position)
#                     print("start_switch",start_switch)

                    # print(gap_counter, base_counter)
                    # print(5,0.2*sequence_length)
#                     print((len(records.seq[:current_position]) - records.seq[:current_position].count('-')), 0.1*sequence_length)
                    if start_switch == 1 and records.seq[current_position] == "-" and gap_counter > 3 and (len(records.seq[:current_position]) - records.seq[:current_position].count('-')) < 0.1*sequence_length:
                        # print("\n\nhere\n\n")
                        start_switch = 0
                        gap_counter = 0
                    
                    if "-" not in records.seq[current_position] and start_switch == 0 :
                        
                        fasta_start_position = current_position
                        start_switch = 1




                    if end_switch == 1 and "-" not in records.seq[current_position] :
                        end_switch = 0
                        gap_counter = 0

#                     print((len(records.seq[current_position:]) - records.seq[current_position:].count('-')), 0.1*sequence_length) 
#                     print((start_switch == 1),records.seq[current_position],(len(records.seq[current_position:]) - records.seq[current_position:].count('-')) < 0.1*sequence_length)
                    if start_switch == 1 and records.seq[current_position] == "-" and end_switch == 0 and (len(records.seq[current_position:]) - records.seq[current_position:].count('-')) < 0.1*sequence_length:
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
            if end_switch == 0:
                end_switch = 1
                fasta_end_position = current_position
            start_switch = 1
    #         break
        print("here here", start_switch,end_switch,min_score_sequence,records.id,fasta_start_position, fasta_end_position )
        if start_switch == 1 and end_switch == 1 and min_score_sequence in records.id :
            alignment_name = records.id
#             print(fasta_start_position, fasta_end_position)
            gene_sequence = records.seq[fasta_start_position:fasta_end_position]
            print(f"{records.id}\n{gene_sequence}")
#     print("_".join(alignment_name.split("_")[:2]).replace("set", "sequence"))
#     clustal_alignment_file = f'alignment_clustal_Error_{error_exon}_translated_genomic_{("_".join(alignment_name.split("_")[:2]).replace("set", "sequence_"))}.fa.txt'
    print(clustal_alignment_file)
    with open(f"{mafft_run_location}{clustal_alignment_file}", 'r') as clustal_file_open:
        return(gene_sequence, alignment_name,clustal_file_open.readlines() )

def run_blast(gene_sequence,query_sequence,alignment_name,output_location,species, error_exon):
    error_exon_location = f"{output_location}/1.Blast_result/{species}/Error_exon_processing/{error_exon}"
    list_of_folders = os.listdir(error_exon_location)
    if "Run_Blast" not in list_of_folders:
        os.mkdir(f"{error_exon_location}/Run_Blast")
    
    with open(f"{error_exon_location}/Run_Blast/local_db.txt",'w') as db_file:
        output = f">genome_fragment\n{gene_sequence}"
        db_file.write(output)
        
    with open(f"{error_exon_location}/Run_Blast/query.txt",'w') as query_file:
        output = f">{alignment_name}\n{query_sequence}"
        query_file.write(output)
    
    makeblast_command = f'cd "{error_exon_location}/Run_Blast/"\nmakeblastdb -in local_db.txt -dbtype nucl'
    subprocess.run(f'{makeblast_command}', shell = True, stderr = subprocess.DEVNULL)
    
    blast_command = f'cd "{error_exon_location}/Run_Blast/"\ntblastn -seg no -query query.txt -db local_db.txt -num_alignments 3 -out blast_out_genome_fragment.htm -html\ntblastn -seg no -query query.txt -db local_db.txt -num_alignments 3 -out blast_out_genome_fragment.txt'
    subprocess.run(f'{blast_command}', shell = True, stderr = subprocess.DEVNULL)
    

def process_genome_fragment_blast_file(output_location,left_oh,right_oh,species,genome_location):
    error_exon_location = f"{output_location}/1.Blast_result/{species}/Error_exon_processing/{error_exon}"
    blast_location = f"{error_exon_location}/Run_Blast"
    
    with open(f"{blast_location}/query.txt", 'r') as query_file:
        query_name_list = [query_file.readlines()[0][1:].rstrip()]
        seq_modi = [[int(left_oh),int(right_oh)]]
        print(seq_modi)

    header = "Species," + "Scaffold," + "Start," + "Stop," + "Complement," + "Error," + "Gene,"+ "Query_start," + "Query_stop,"+ "Query_Length\n" 
    Output_Sequence = header
    scaff = "Intial_value"
    scaff_old = "Intial_value"
    old_end = 0
#     species_name = annotated_species_name

    for i in range(len(query_name_list)):
        query_name = query_name_list[i]
        print(query_name)
#         return
        Length_switch = "0"
        
        with open(f"{blast_location}/blast_out_genome_fragment.txt",'r') as tblast_out:
            lines_in_file = tblast_out.readlines()

        result_section_switch = 0
        start_coor_switch = 0
        query_start_coor_switch = 0
        stop_coor_switch = 0
        error = "N"
        break_switch = 0

        start = 0
        stop = 0
        start_coor = 0
        stop_coor = 0
        query_length = 0
        gt_ag = "N"

        for lines in lines_in_file:

#             print(lines)
            if query_name in lines:
            #Initialize that results can now be checked
                result_section_switch = 1
                query_species_split = lines.split(" ")[1].split("_")
                query_species = str(query_species_split[1]+"_"+query_species_split[2].rstrip())

            if result_section_switch == 1 and "Lambda" in lines:
            #This block indicates end of the results block in blast output
                result_section_switch == 0
                
                break

            if result_section_switch == 1:
            #While checking the result
                if "Length=" in lines and Length_switch == "0":
                #Get query length from the blast output
                    
                    query_length = int(lines.split("=")[1].rstrip())
                    
                    Length_switch = 1 #Indicated length has been acquired
                    
                if ("Score" in lines or ">" in lines) and (start_coor_switch == 1):
    #                print (lines)
                    break
        
                if ">" in lines:
                #Start of the first result
                    scaff = lines.split(" ")[1] #Scaffold from the result
                    if scaff_old != "Intial_value" and scaff_old != scaff:
                        error = "Y"
                    scaff_old = scaff
                    
                if "Query" in lines and "=" not in lines:
                #Read the query line in output
                    if query_start_coor_switch == 0:
#                        print(lines)
                        query_start_coor = int(lines.split(" ")[2])
                        query_start_coor_switch = 1
                        #Query start coordinate fixed
            
                    query_stop_coor =int(lines.split(" ")[-1][:-1])
                    #Keep getting query stop coordinates for multiline result
    #                print (stop_coor)
                    
                if "Sbjct" in lines:
                #Read the blast target line
                    if start_coor_switch == 0:
                        start_coor = int(lines.split(" ")[2])
                        start_coor_switch = 1
                    stop_coor =int(lines.split(" ")[-1][:-1])
                    #Keep getting target stop coordinates for multiline result
                
                

        
        if break_switch == 1:
            break
        print(start_coor,stop_coor)
        if start_coor < stop_coor:
            complement = "0" #Forward complement
            
            length = (stop_coor-start_coor)/3
            start = start_coor
            stop = stop_coor

        elif start_coor > stop_coor:
            complement = "1" #Reverse complement
            length = (-stop_coor+start_coor)/3
            start = stop_coor
            stop = start_coor

        else:
            error = "Y"
        
        
        
        
        seq_length = query_length
        if (start != 0 or stop != 0):
            start_modifier = seq_modi[i][0]
            stop_modifier = seq_modi[i][1]
        else:
            start_modifier = 0
            stop_modifier = 0  
        
        
#Adding or removing 3' and 5' overhangs for forward and reverse complement
    #For forward complement
        if complement == "0":
            start = int(start) - int(start_modifier)
            stop = int(stop) +  int(stop_modifier)
            if old_end != 0 and old_end > stop:

                error = "Y"
            old_end = stop

    #For reverse complement
        if complement == "1":
            start = int(start) - int(stop_modifier)
            stop = int(stop) +  int(start_modifier)
            if old_end != 0 and old_end < stop:
                error = "Y"
            old_end = stop

    #Simple check for lenghth
        if start == 0 or stop == 0:
            error = "Y"

        genome_file = SeqIO.parse(f"{blast_location}/local_db.txt", 'fasta')
        print("reached here")
        
        for records in genome_file:
            old_start = start
            old_stop = stop
            ag = "N"
            gt = "N"
            stop_counter = 0
            while True:
                print(f"sequence:\n{records.seq[start+start_modifier-1:stop]}")
                translated_sequence = records.seq[start+start_modifier-1:stop].translate()
                print(f"sequence:\n{translated_sequence}")
                if "*" in translated_sequence:
                    stop_counter +=1
                    if ag == "N":
                        start = old_start + 3*stop_counter
                    if gt =="N":
                        stop = old_stop - 3*stop_counter
                if stop_counter > 90:
                    return ("Error!! Too many stops")
                print(f"left = {records.seq[start-3:start-1]}, right ={(records.seq[stop:stop + 2])}"  ), 
                if (records.seq[start-3:start-1]).lower() == "ag" and ag != "Y":
                    
                    ag = "Y"
                    
                elif ag != "Y":
                    start -= 3
                
                if (records.seq[stop:stop + 2]).lower() == "gt" and gt != "Y":
                    gt = "Y"
                elif gt != "Y":
                    stop +=3
                if old_start - start > 1000 or stop - old_stop > 1000:
                    break
                if gt == "Y" and ag == "Y":
                    gt_ag = "Y"
                    break
        
        query_location = f"{blast_location}/new_query_spliced.txt"
        with open(query_location , 'w') as query_file_new:

            sequence_translated = records.seq[start+start_modifier-1:stop].translate()
            print(sequence_translated)
            # proceed_test = input("Proceed with this?")
            # while True:
            #     if proceed_test.lower()[0] == "n":
            #         assert False
            #     elif proceed_test.lower()[0] == "y":
            #         break
            if "*" in sequence_translated:
                print("Errror in Spliced query")
                assert False
            
            output = f">{query_name}\n{sequence_translated}"
            query_file_new.write(output)

        genome_file = get_genome_file(genome_location,species)
        genome = f"{genome_location}/{species}/{genome_file}"
        out_location = f"{blast_location}"
    
        
        blast_command = f'tblastn -seg no -query "{query_location}" -db "{genome}" -num_alignments 3 -out "{out_location}/blast_out.htm" -html'
        # print(blast_command)
        # subprocess.run(f'{blast_command}', shell = True, stderr = subprocess.DEVNULL)
        subprocess.run(f'{blast_command}', shell = True, stderr = subprocess.DEVNULL)
    
        blast_command = f'tblastn -seg no -query "{query_location}" -db "{genome}" -num_alignments 3 -out "{out_location}/blast_out.txt"'
        subprocess.run(f'{blast_command}', shell = True, stderr = subprocess.DEVNULL)
    

                
        # acceptor, donor, don_line, acc_line = process_spiceator_result(start, stop,annotated_genome_location,annotated_species_name)

#         if acceptor == "Y" and donor  == "Y":
#             splice_prediction = "Y"
#         else:
#             splice_prediction = "N"
            
        return(gt_ag) 
    
    

def process_genome_blast_file(output_location,left_oh,right_oh,species,genome_location,ag_gt,query_name_original):
    error_exon_location = f"{output_location}/1.Blast_result/{species}/Error_exon_processing/{error_exon}"
    blast_location = f"{error_exon_location}/Run_Blast"

    with open(f"{blast_location}/new_query_spliced.txt", 'r') as query_file:
        query_name_list = [query_file.readlines()[0][1:].rstrip()]
        seq_modi = [[int(left_oh),int(right_oh)]]
        print(seq_modi)

    header = "Species," + "Scaffold," + "Start," + "Stop," + "Complement," + "Error," + "Gene,"+ "Query_start," + "Query_stop,"+ "Query_Length," +  "AG_GT," + "Spliceator_prediction\n"
    
    Output_Sequence = header
    scaff = "Intial_value"
    scaff_old = "Intial_value"
    old_end = 0
    species_name = species

    for i in range(len(query_name_list)):
        query_name = query_name_list[i]
        Length_switch = "0"
        
        with open(f"{blast_location}/blast_out.txt",'r') as tblast_out:
            lines_in_file = tblast_out.readlines()

        result_section_switch = 0
        start_coor_switch = 0
        query_start_coor_switch = 0
        stop_coor_switch = 0
        error = "N"
        break_switch = 0

        start = 0
        stop = 0
        start_coor = 0
        stop_coor = 0
        query_length = 0
        gt_ag = "N"

        for lines in lines_in_file:

            # print(lines)
            # print(query_name)
            if query_name in lines:
            #Initialize that results can now be checked
                result_section_switch = 1
                query_species_split = lines.split(" ")[1].split("_")
                query_species = str(query_species_split[1]+"_"+query_species_split[2].rstrip())

            if result_section_switch == 1 and "Lambda" in lines:
            #This block indicates end of the results block in blast output
                result_section_switch == 0
                
                break

            if result_section_switch == 1:
            #While checking the result
                if "Length=" in lines and Length_switch == "0":
                #Get query length from the blast output
                    
                    query_length = int(lines.split("=")[1].rstrip())
                    
                    Length_switch = 1 #Indicated length has been acquired
                    
                if ("Score" in lines or ">" in lines) and (start_coor_switch == 1):
    #                print (lines)
                    break
        
                if ">" in lines:
                #Start of the first result
                    scaff = lines.split(" ")[1] #Scaffold from the result
                    if scaff_old != "Intial_value" and scaff_old != scaff:
                        error = "Y"
                    scaff_old = scaff
                    
                if "Query" in lines and "=" not in lines:
                #Read the query line in output
                    if query_start_coor_switch == 0:
#                        print(lines)
                        query_start_coor = int(lines.split(" ")[2])
                        query_start_coor_switch = 1
                        #Query start coordinate fixed
            
                    query_stop_coor =int(lines.split(" ")[-1][:-1])
                    #Keep getting query stop coordinates for multiline result
    #                print (stop_coor)
                    
                if "Sbjct" in lines:
                #Read the blast target line
                    if start_coor_switch == 0:
                        start_coor = int(lines.split(" ")[2])
                        start_coor_switch = 1
                    stop_coor =int(lines.split(" ")[-1][:-1])
                    #Keep getting target stop coordinates for multiline result
                
                

        
        if break_switch == 1:
            break
        
        print(f"start_coordinate : {start_coor},stop_coordinate : {stop_coor}")
        if start_coor < stop_coor:
            complement = "0" #Forward complement
            
            length = (stop_coor-start_coor)/3
            start = start_coor
            stop = stop_coor

        elif start_coor > stop_coor:
            complement = "1" #Reverse complement
            length = (-stop_coor+start_coor)/3
            start = stop_coor
            stop = start_coor

        else:
            error = "Y"
        
        
        
        
        seq_length = query_length
        if (start != 0 or stop != 0):
            start_modifier = seq_modi[i][0]
            stop_modifier = seq_modi[i][1]
        else:
            start_modifier = 0
            stop_modifier = 0  
        #Check if the length of target (blast hit) is significantly smaller than query
        
#Adding or removing 3' and 5' overhangs for forward and reverse complement
    #For forward complement
        if complement == "0":
            start = int(start) - int(start_modifier)
            stop = int(stop) +  int(stop_modifier)
            if old_end != 0 and old_end > stop:

                error = "Y"
            old_end = stop

    #For reverse complement
        if complement == "1":
            start = int(start) - int(stop_modifier)
            stop = int(stop) +  int(start_modifier)
            if old_end != 0 and old_end < stop:
                error = "Y"
            old_end = stop

    #Simple check for lenghth
        if start == 0 or stop == 0:
            error = "Y"

        splice_prediction = "Y"   
        output_format = str(species_name)+"," + str(scaff) +"," + str(start)+"," + str(stop)+"," + str(complement)+"," + str(error)+  ","+ str(query_name_original)+","+ str(query_start_coor)+","+str(query_stop_coor)+","+str(query_length)+ "," + ag_gt + "," + splice_prediction +"\n"  
        # print(Output_Sequence)
        return(output_format)
    
    

with open(f"{output_location}/error_exons.txt",'r') as error_file:
    error_file_lines = error_file.readlines()
    
for line in error_file_lines:
    if len(line) != 1:
        print(line)
        line_split = line.split(",")
        species = line_split[0]
        error_exon = line_split[6]
        
        make_folder_error_exon(output_location,species,error_exon)
        left_oh,right_oh = make_query(output_location,error_exon,species)
        fragment_start, fragment_end, scaffold, complement,query_name = get_genomic_coordinates(output_location,error_exon,species)
        print(f"Gene_start = {fragment_start}, Gene_end = {fragment_end}")
        gene_sequence = get_gene_sequence(genome_location, 
                          species,                       
                          scaffold,
                          fragment_start,
                          fragment_end,
                          complement,
                          output_location)
        
        make_raw_files_for_alignment(gene_sequence,output_location,species,error_exon)
        mafft_run_location = run_mafft(output_location,species,error_exon)
        
        query_sequence, alignment_name,alignment_file = process_mafft_output(mafft_run_location, error_exon)
        
        print(f'Alignment File: {"".join(alignment_file)}')
#        query_check = input("Query OK?")
        query_check = "y"
        if query_check[0].lower() == 'n':
            query_sequence = input("Add new query :")
            if query_sequence == '':
                print(species, error_exon)
                assert False
        run_blast(gene_sequence,query_sequence,alignment_name,output_location,species, error_exon)
        ag_gt = process_genome_fragment_blast_file(output_location,left_oh,right_oh,species,genome_location)
        if ag_gt != "Error!! Too many stops":
            coordinate_output = process_genome_blast_file(output_location,left_oh,right_oh,species,genome_location,ag_gt,query_name)
            print(coordinate_output)
        else:
            print(ag_gt)
        print("\n\n\nFix Overhang!! Proceed?\n\n\n")
#         print(gene_sequence)
        
        