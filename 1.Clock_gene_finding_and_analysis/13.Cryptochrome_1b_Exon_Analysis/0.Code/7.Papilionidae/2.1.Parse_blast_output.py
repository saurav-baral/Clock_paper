# -*- coding: utf-8 -*-
"""
Created on Tue May 28 13:35:56 2024

@author: sauba
"""
import os
from Bio import SeqIO

from Bio.Seq import Seq


def sequence_extractor(Species, Scaff, reverse_c, start, end, frame = 1, trans = 1 ):
    from Bio import SeqIO
    from Bio.Seq import Seq
    import os
    
    out_Seq = ''
    
    
    
    #Genome folder
    # entries = os.listdir("/mnt/f/Genomes_2023-12-26/"+Species)
    entries = os.listdir("/mnt/griffin/saubar/Genomes_2023-12-26/"+Species)
    
    #Get genome file from Genome folder
    for file_names  in entries:
        if ".nhr" in file_names:
            Genome_name = file_names[:-4]
            break
    
    #Read the genome file
    # fasta_file = SeqIO.parse(("/mnt/f/Genomes_2023-12-26/"+Species+"/"+Genome_name),'fasta')
    fasta_file = SeqIO.parse(("/mnt/griffin/saubar/Genomes_2023-12-26/"+Species+"/"+Genome_name),'fasta')
    
    #Extract the sequence 
    for record in fasta_file:
        if record.id == Scaff:
            sequence = str(record.seq)
            out_Seq = Seq(sequence[start-1:end])
            if reverse_c == 1:
                out_Seq = out_Seq.reverse_complement()
            break
                
    if len(out_Seq) < 10000: #fixing error due to mistake in typing the coordinate, change this for longer sequence
#        print (out_Seq)
        if frame == 1 :
            out_trans = out_Seq[0:]
        if frame == 2 :
            out_trans = out_Seq[1:]
        if frame == 3 :
            out_trans = out_Seq[2:]
    else:
        print ("too long")
        assert(False)
        
        
    #    break
    
    if trans == 1:
        return (out_trans.translate())
    else:
        return(out_Seq)

def get_sequence(Species, Scaff, reverse_c, start, end, frame = 1, trans = 1 ):
    
    #Genome folder
    entries = os.listdir("/mnt/griffin/saubar/Genomes_2023-12-26/"+Species)
    
    #Get genome file from Genome folder
    for file_names  in entries:
        if ".nhr" in file_names:
            Genome_name = file_names[:-4]
            break
    
    #Read the genome file
    fasta_file = SeqIO.parse(("/mnt/griffin/saubar/Genomes_2023-12-26/"+Species+"/"+Genome_name),'fasta')
    
    #Extract the sequence 
    for record in fasta_file:
        if record.id == Scaff:
            sequence = str(record.seq)
#             out_Seq = Seq(sequence[start-1:end])
#             if reverse_c == 1:
#                 sequence = out_Seq.reverse_complement()
            return(sequence)
    
#


file_location = os.getcwd()
family_group = file_location.split("/")[-1]

blast_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/13.Cryptochrome_1b_Exon_Analysis/{family_group}"
query_location = f"{blast_location}/1.Query"
query_species = "Papilio_xuthus"
query_transcript_list = os.listdir(f"{query_location}/{query_species}")
if "desktop.ini" in query_transcript_list:
    query_transcript_list.remove("desktop.ini")

query_transcript = query_transcript_list[0]
#Get Species list, the list may be a list of files or a single species

species_list = os.listdir(f"{blast_location}/1.Blast_result")
# species_list = ["Pararge_aegeria"]
if "desktop.ini" in species_list:
    species_list.remove('desktop.ini')

#get list of query names

query_file = SeqIO.parse(f"{blast_location}/1.Query/{query_species}/{query_transcript}/query_protein.fa", 'fasta')
query_name_list = []

for records in query_file:
    query_name_list.append(records.id)



for species in species_list:
    # print(species_name)
    # species = species_name.split(".")[1]
    header = "Species," + "Scaffold," + "Start," + "Stop," + "Complement," + "Error," + "Gene,"+ "Query_start," + "Query_stop,"+ "Query_Length\n" 
    Output_Sequence = header
    scaff = "Intial_value"
    scaff_old = "Intial_value"
    old_end = 0
    
#Run for each query in the query list
    for i in range(len(query_name_list)):
        query_name = query_name_list[i]
#         print(f"Query name = {query_name}")
        Length_switch = "0"
        with open(f"{blast_location}/1.Blast_result/{species}/{species}_blast_out.txt",'r') as tblast_out:
            lines_in_file = tblast_out.readlines()
        #print(lines_in_file)

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

        for lines in lines_in_file:

#             print(lines)
            if query_name in lines:
#                 print(lines)
#                 assert False
            #Initialize that results can now be checked
                result_section_switch = 1
#                 query_species_split = lines.split(" ")[1].split("_")
#                 query_species = str(query_species_split[0].split(".")[1]+"_"+query_species_split[1].rstrip())
                print(query_species)
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
            length = 0
            start = 000
            stop = 000
            query_start_coor = 000
            query_stop_coor = 000
            error = "Y"
            complement  = 0
        
        
        
        

        #Check if the length of target (blast hit) is significantly smaller than query
        if length < query_length - 0.2*query_length:
            error = "Y"

        old_trans = ''
#             print(start, stop)
        
    
    
        

#Check for Met or STOP at the start and end
        sequence = ''
        while True and (start != 0 or stop != 0):
            if (query_name == query_name_list[0] or query_name == query_name_list[-1]):
                frame = 1
#                 translated_sequence = (sequence_extractor(species, scaff, int(complement), start, stop,frame ))
                if sequence == '':
                    sequence = get_sequence(species, scaff, int(complement), start, stop,frame)
                if int(complement) == 1:
                    translated_sequence = Seq(sequence[start-1:stop]).reverse_complement().translate()
                else:
                    translated_sequence = Seq(sequence[start-1:stop]).translate()
                if translated_sequence == old_trans:
                    break
                    
#Check for Met at the beginning of first exon

                if query_name == query_name_list[0] :
                    if translated_sequence[0] != "M":
                        print(translated_sequence)
#                             print(start, stop)
                        if complement == "0":
                            start = int(start) - 3
                        if complement == "1":
                            stop = int(stop) + 3
                        old_trans = translated_sequence
                    if "*" in translated_sequence:
                        error = "Y"
                        break
                    if translated_sequence[0] == "M":
                        break

#Check for stop at the end of last exon

                if query_name == query_name_list[-1]:
                    if translated_sequence[-1] != "*":
                        print(translated_sequence)
#                             print(start, stop)
                        if complement == "0":
                            stop = int(stop) + 3
                        if complement == "1":
                            start = int(start) - 3
                        old_trans = translated_sequence
                    if translated_sequence[-1] == "*":

                        break
            else:
                break
        
        
#Add offset at the beginning and end of each exon
        if (start != 0 or stop != 0):
            
            start_modifier = int(query_name.split("Frame")[1][1])
            stop_modifier = int(query_name.split("rightoh")[1][1])
        else:
            start_modifier = 0
            stop_modifier = 0            

        seq_length = query_length
        
        
#Get exon start and end by adding or removing codons
        
    #if not the beginning
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
            
           
        output_format = str(species.split("\n")[0])+"," + str(scaff) +"," + str(start)+"," + str(stop)+"," + str(complement)+"," + str(error)+  ","+ str(query_name)+","+ str(query_start_coor)+","+str(query_stop_coor)+","+str(query_length)+ "\n"  
        print(output_format)
#         assert False
        Output_Sequence = Output_Sequence + output_format

    print(Output_Sequence)
   
    output_file = open(f"{blast_location}/1.Blast_result/{species}/{species}_coordinates_old.csv",'w')
    output_file.write(Output_Sequence)
    output_file.close()

