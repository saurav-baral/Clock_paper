# -*- coding: utf-8 -*-
"""
Created on Fri May 31 17:03:19 2024

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


annotated_genome_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/15.Timeout_Exon_Analysis/{family_group}/1.Blast_result"
query_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/15.Timeout_Exon_Analysis/{family_group}/1.Query"
genome_location = "/mnt/griffin/saubar/Genomes_2023-12-26"

species_list = [f"{sys.argv[1]}"]
if "desktop.ini" in species_list:
    species_list.remove("desktop.ini")


def run_blast_with_new_query(annotated_genome_location,
                             annotated_species_name,
                             error_exon, 
                            query_species,
                             genome_location, 
                             species):
    
    genome_file = get_genome_file(genome_location,species)
    query_location = f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_blast/new_query.txt"
    genome = f"{genome_location}/{species}/{genome_file}"
    out_location = f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_blast"
    try:
        with open(f"{out_location}/new_query.txt", 'r') as query_file_open:
            query_file_lines = query_file_open.readlines()
    except:
        return ("Query_error")
    if (len(query_file_lines)) == 1:
        return ("Query_error")
    print(len(query_file_lines[1])) 
    if (len(query_file_lines[1])) < 5:
        return ("Query_error")
    # print(len(query_file_lines[1]))
    # genome_fragment_out = f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_blast/local_genomic_fragment.fa"
    # with io.open(genome_fragment_out,'w') as out_file:
    #     output = f">Genome_fragment_{error_exon}\n{gene_sequence}"
    #     out_file.write(output)

    local_genomic_fragment_location = f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon/"
    # cd_command = f'cd "{local_genomic_fragment_location}"\nmakeblastdb -in gene_sequence_all.fa -dbtype nucl\n'
    # # os.system(f'{cd_command}')
    # subprocess.run(f'{cd_command}', shell = True, stderr = subprocess.DEVNULL)
    # # print(mkdb_command)
    # # os.system(f'{mkdb_command}')
    blast_command = f'cd "{out_location}"\ntblastn -seg no -query new_query.txt -db ../../../Period_gene_genomic_sequence_individual_exon/gene_sequence_all.fa -num_alignments 3 -out blast_out_genome_fragment.htm -html'
    # print(blast_command)
    # os.system(f'{blast_command}')
    subprocess.run(f'{blast_command}', shell = True, stderr = subprocess.DEVNULL)
    # subprocess.run(f'{blast_command}', shell = True)
    blast_command = f'cd "{out_location}"\ntblastn -seg no -query new_query.txt -db ../../../Period_gene_genomic_sequence_individual_exon/gene_sequence_all.fa -num_alignments 3 -out blast_out_genome_fragment.txt'
    # os.system(f'{blast_command}')
    subprocess.run(f'{blast_command}', shell = True, stderr = subprocess.DEVNULL)


def get_genome_file(genome_location,species):
    list_of_files_in_genome_folder = os.listdir(f"{genome_location}/{species}")
    for file in list_of_files_in_genome_folder:
        if file.endswith("_genomic.fna"):
            genome_file = file
    return(genome_file)  


def get_query_info(query_location,query_species,query_transcript,query_exon):
    with open(f"{query_location}/{query_species}/{query_transcript}/query_{query_exon}.fa", 'r') as query_file:
        query_file_list = query_file.readlines()
        query_fasta_sequence = "".join(query_file_list)
        left_overhang = query_file_list[0].split("Frame")[1][1]
        right_overhang = query_file_list[0].split("rightoh")[1][1]
        original_query_name = query_file_list[0]
    return(left_overhang,right_overhang,original_query_name)


def process_genome_fragment_blast_file(annotated_genome_location, 
                                       annotated_species_name,
                                       error_exon,
                                       left_overhang,
                                       right_overhang,                                                                                               
                                       original_query_name,
                                      query_species_original ):
    blast_location = f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species_original}/{error_exon}/for_blast"
    
    with open(f"{blast_location}/new_query.txt", 'r') as query_file:
        query_name_list = [query_file.readlines()[0][1:].rstrip()]
        seq_modi = [[int(left_overhang),int(right_overhang)]]
        print(seq_modi)

    header = "Species," + "Scaffold," + "Start," + "Stop," + "Complement," + "Error," + "Gene,"+ "Query_start," + "Query_stop,"+ "Query_Length\n" 
    Output_Sequence = header
    scaff = "Intial_value"
    scaff_old = "Intial_value"
    old_end = 0
    species_name = annotated_species_name

    for i in range(len(query_name_list)):
        query_name = query_name_list[i]
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
        if length < query_length - 0.2*query_length:
            error = "Y"

        old_trans = ''

        if query_start_coor != "1" and query_name != query_name_list[0]:
            if complement == "0":
                start = int(start) - 3*(int(query_start_coor)-1)                
            if complement == "1":
                stop = int(stop) + 3*(int(query_start_coor)-1)
    
    #For the end
        if query_stop_coor != str(seq_length) and query_name != query_name_list[-1]:
            if complement == "0":
                stop = int(stop) + 3*(int(seq_length)-int(query_stop_coor))
            if complement == "1":
                

                start = int(start) - 3*(int(seq_length)-int(query_stop_coor))
 
        
        
        
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

        genome_file = SeqIO.parse(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon/gene_sequence_all.fa", 'fasta')
        print("reached here")
        print(f"Before splice: {start},{stop}")
        acceptor, donor, don_line, acc_line = process_spiceator_result(start, stop,annotated_genome_location,annotated_species_name,query_species_original)
        print(f"after splice: {start},{stop}")
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
                print(f"left = {records.seq[start-3:start-1]}, right ={(records.seq[stop:stop + 2])}, stop_counter = {stop_counter}"  ), 
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
                if stop_counter > 10:
                    gt = "Y"
                    ag = "Y"
                
        
        query_location = f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species_original}/{error_exon}/for_blast/new_query_spliced.txt"
        with open(query_location , 'w') as query_file_new:
            # sequence_translated = records.seq[start+start_modifier-1:stop]
            # print(sequence_translated)
            
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
            
            output = f">{original_query_name[1:]}\n{sequence_translated}"
            query_file_new.write(output)

        genome_file = get_genome_file(genome_location,species)
        genome = f"{genome_location}/{species}/{genome_file}"
        out_location = f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species_original}/{error_exon}/for_blast"
    
        
        blast_command = f'tblastn -seg no -query "{query_location}" -db "{genome}" -num_alignments 3 -out "{out_location}/blast_out.htm" -html'
        # print(blast_command)
        # subprocess.run(f'{blast_command}', shell = True, stderr = subprocess.DEVNULL)
        subprocess.run(f'{blast_command}', shell = True, stderr = subprocess.DEVNULL)
    
        blast_command = f'tblastn -seg no -query "{query_location}" -db "{genome}" -num_alignments 3 -out "{out_location}/blast_out.txt"'
        subprocess.run(f'{blast_command}', shell = True, stderr = subprocess.DEVNULL)
    

                
        # acceptor, donor, don_line, acc_line = process_spiceator_result(start, stop,annotated_genome_location,annotated_species_name)

        if acceptor == "Y" and donor  == "Y":
            splice_prediction = "Y"
        else:
            splice_prediction = "N"
            
        return(start,stop,gt_ag, splice_prediction) 
        # output_format = str(species_name.split("\n")[0])+"," + str(scaffold) +"," + str(start)+"," + str(stop)+"," + str(complement)+"," + str(error)+  ","+ str(query_name)+","+ str(query_start_coor)+","+str(query_stop_coor)+","+str(query_length)+ "\n"  

def process_spiceator_result(start_coordinate, stop_coordinate,annotated_genome_location,annotated_species_name,query_species_original):
    results_location = f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species_original}/"
#     list_of_files_here = os.listdir(results_location)
#     spliceator_results_file = ''
#     for files in list_of_files_here:
#         if files.startswith("Spliceator_results"):
#             spliceator_results_file = files
#     if spliceator_results_file == '':
#         print("Splice file missing")
#         assert False

#     with open(f"{results_location}/{spliceator_results_file}", 'r') as splice_file:
#         splice_file_content = splice_file.readlines()

    splice_acceptor_presence = "N"
    splice_donor_presence = "N"
    splice_donor = ''
    splice_acceptor = ''
    start_coordinate = 0
    stop_coordinate = 0
#     for lines in splice_file_content:
#         line_split = lines.split(";")
#         if line_split[0]=="Acceptor" and splice_acceptor_presence == "N":
#             acceptor_start = int(line_split[1])
#             acceptor_end = acceptor_start + len(line_split[3])
            
#             if start_coordinate >=acceptor_start and start_coordinate <=acceptor_end:
#                 splice_acceptor_presence = "Y"
#                 # while True:
#                 #     if start_coordinate >= acceptor_start:
#                 #         start_coordinate -= 3
#                 #     else:
#                 #         break
#                 splice_acceptor = lines
#         if line_split[0]=="Donor" and splice_donor_presence == "N":
#             donor_start = int(line_split[1])
#             donor_end = donor_start + len(line_split[3])
#             if stop_coordinate >=donor_start and stop_coordinate <=donor_end:
#                 splice_donor_presence = "Y"
#                 # while True:
#                 #     if stop_coordinate >=donor_start:
#                 #         stop_coordinate -= 3
#                 #     else:
#                 #         break
#                 splice_donor = lines

    return(splice_acceptor_presence, splice_donor_presence, splice_donor, splice_acceptor)
                
def process_genome_blast_file(annotated_genome_location, annotated_species_name,error_exon,left_overhang,right_overhang, ag_gt, splice_prediction,original_query_name, query_species ):
    blast_location = f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}/{error_exon}/for_blast"
    with open(f"{blast_location}/new_query_spliced.txt", 'r') as query_file:
        query_name_list = [query_file.readlines()[0][1:].rstrip()]
        seq_modi = [[int(left_overhang),int(right_overhang)]]
        print(seq_modi)

    header = "Species," + "Scaffold," + "Start," + "Stop," + "Complement," + "Error," + "Gene,"+ "Query_start," + "Query_stop,"+ "Query_Length," +  "AG_GT," + "Spliceator_prediction\n"
    
    Output_Sequence = header
    scaff = "Intial_value"
    scaff_old = "Intial_value"
    old_end = 0
    species_name = annotated_species_name

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
        if length < query_length - 0.2*query_length:
            error = "Y"

        old_trans = ''

        if query_start_coor != "1" and query_name != query_name_list[0]:
            if complement == "0":
                start = int(start) - 3*(int(query_start_coor)-1)                
            if complement == "1":
                stop = int(stop) + 3*(int(query_start_coor)-1)
    
    #For the end
        if query_stop_coor != str(seq_length) and query_name != query_name_list[-1]:
            if complement == "0":
                stop = int(stop) + 3*(int(seq_length)-int(query_stop_coor))
            if complement == "1":
                

                start = int(start) - 3*(int(seq_length)-int(query_stop_coor))
 
        
        
        
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

            
        output_format = str(species_name)+"," + str(scaff) +"," + str(start)+"," + str(stop)+"," + str(complement)+"," + str(error)+  ","+ str(query_name)+","+ str(query_start_coor)+","+str(query_stop_coor)+","+str(query_length)+ "," + ag_gt + "," + splice_prediction +"\n"  
        # print(Output_Sequence)
        return(output_format)





blast_error = []
for annotated_species_name in species_list:
    # annotated_species_name = "Lasiommata_megera"
    species = annotated_species_name
    

    list_of_query_species = os.listdir(query_location)
    if "desktop.ini" in list_of_query_species:
        list_of_query_species.remove("desktop.ini")

    for query_species in list_of_query_species:
    
        query_transcript_list = os.listdir(f"{query_location}/{query_species}")
        if "desktop.ini" in query_transcript_list:
            query_transcript_list.remove("desktop.ini")

        query_transcript = query_transcript_list[0]

        

        list_of_exons_folders = os.listdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon_{query_species}")

        error_exon_list = []
        for folders in list_of_exons_folders:

            if folders.startswith("Exon"):
                if int(folders.split("_")[1]) > 1:
                    error_exon_list.append(folders)
        print(error_exon_list)
        # break


    #     error_exon_list = ["Exon_5"]
        for error_exon in error_exon_list:
            print(annotated_species_name)
            print(f"Processing {error_exon}")
            coordinate_location = f"{annotated_genome_location}/{species}"
            files_in_coordinate_location = os.listdir(coordinate_location)
            coordinate_backup = 0
            query_coordinate = 0
            # print(files_in_coordinate_location)
            for file_names in files_in_coordinate_location:
                if file_names.endswith("_coordinates_old.csv"):
                    coordinate_backup = 1
                if file_names.endswith(f"_coordinates_{query_species}.csv"):
                    print(file_names)
                    query_coordinate = 1
            print(query_coordinate, "query coordinate")
            if query_coordinate == 0:
                print(file_names)
                copy_command = f'cp "{coordinate_location}/{species}_coordinates_old.csv" "{coordinate_location}/{species}_coordinates_{query_species}.csv"'
                subprocess.run(f'{copy_command}', shell = True, stderr = subprocess.DEVNULL)
            if coordinate_backup == 0:
                copy_command = f'cp "{coordinate_location}/{species}_coordinates.csv" "{coordinate_location}/{species}_coordinates_old.csv"'
                subprocess.run(f'{copy_command}', shell = True, stderr = subprocess.DEVNULL)

            query_error = run_blast_with_new_query(annotated_genome_location,
                                         annotated_species_name,
                                         error_exon, 
                                         query_species,
                                         genome_location, 
                                         species )

            if query_error == "Query_error":
                print(query_error)
                blast_error.append(f"{annotated_species_name},{error_exon}")
                continue

            left_overhang,right_overhang,original_query_name = get_query_info(query_location,query_species,query_transcript,error_exon)

            try:
                start_coordinate,stop_coordinate,gt_ag, splice_prediction = process_genome_fragment_blast_file(annotated_genome_location, 
                                                                                                           annotated_species_name,
                                                                                                           error_exon,
                                                                                                           left_overhang,
                                                                                                           right_overhang,                                                                                               
                                                                                                           original_query_name,
                                                                                                          query_species)

                new_coordinate_file_line = process_genome_blast_file(annotated_genome_location, 
                                                             annotated_species_name,
                                                             error_exon,
                                                             left_overhang,
                                                             right_overhang,
                                                             gt_ag, 
                                                             splice_prediction,
                                                             original_query_name,
                                                                query_species)
            except:
                print("Error in Blast")
                blast_error.append(f"{annotated_species_name},{error_exon}")
                continue

            print(start_coordinate,stop_coordinate,gt_ag, splice_prediction)
            print(new_coordinate_file_line)
            output = ''
            with open(f"{coordinate_location}/{species}_coordinates_{query_species}.csv", 'r') as open_coordinate_file:
                coordinate_file_list = open_coordinate_file.readlines()
                # print(coordinate_file_list)
                # assert False
            for lines in coordinate_file_list:
                if lines.split(",")[6].endswith(error_exon):
                    lines = new_coordinate_file_line
                output += lines

            with open(f"{coordinate_location}/{species}_coordinates_{query_species}.csv", 'w') as out_coordinate_file:
                out_coordinate_file.write(output)


print("\n".join(blast_error))



