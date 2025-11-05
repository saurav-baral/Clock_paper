# -*- coding: utf-8 -*-
"""
Created on Fri May 31 15:59:52 2024

@author: sauba
"""

import os
import io
import subprocess
import sys



family_group = "3.Satyrine_Exon5_split_for_dup"
genome_location =  "/mnt/griffin/saubar/Genomes_2023-12-26"

blast_output_location = f"/mnt/griffin/saubar/Working_folder_circadian_rhythm/10.Cycle_Exon_Analysis/{family_group}/1.Blast_result"

species_list = [f"{sys.argv[1]}"]

if "desktop.ini" in species_list:
    species_list.remove("desktop.ini")


def get_gene_sequence(genome_location, species, genome_file, scaffold,gene_start,gene_end,complement,annotated_genome_location,annotated_species_name):
    from Bio import SeqIO
   
    list_of_files_inside_annotated_species_folder = os.listdir(f"{annotated_genome_location}/{annotated_species_name}")

    if "Period_gene_genomic_sequence_individual_exon" not in list_of_files_inside_annotated_species_folder:
        os.mkdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon")

    list_of_files_inside_indiv_exon_folder = os.listdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon")
    if "temp" not in list_of_files_inside_indiv_exon_folder:
         os.mkdir(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon/temp")

    os.system(f'samtools faidx "{genome_location}/{species}/{genome_file}"')
    os.system(f'samtools faidx "{genome_location}/{species}/{genome_file}" {scaffold}:{gene_start}-{gene_end} > "{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon/temp/temp_genome.fa"')
    genome = SeqIO.parse(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon/temp/temp_genome.fa", "fasta")
    
    for entries in genome:
        gene_sequence = entries.seq
        if complement == "1":
            gene_sequence = gene_sequence.reverse_complement()
        break
    return (gene_sequence)

def get_genome_file(genome_location,species):
    list_of_files_in_genome_folder = os.listdir(f"{genome_location}/{species}")
    for file in list_of_files_in_genome_folder:
        if file.endswith("_genomic.fna"):
            genome_file = file
    return(genome_file)   

def get_sequence(blast_output_location,genome_location, species, genome_file,annotated_genome_location,annotated_species_name):
    list_of_files_in_species_folder = os.listdir(f"{blast_output_location}/{species}")
    coordinate_file_name = ''
    for file_names in list_of_files_in_species_folder:
        if file_names.endswith("_coordinates_old.csv"):
            coordinate_file_name = file_names

    with io.open(f"{blast_output_location}/{species}/{coordinate_file_name}", 'r') as temp_file_open:
        coordinate_file_lines = temp_file_open.readlines()

    upstream_exon_line =coordinate_file_lines[1].split(",")
    if upstream_exon_line[5] == "Y":
        print("First Exon Error!!")
        assert False

    downstream_exon_line = coordinate_file_lines[-1].split(",")
    if downstream_exon_line[5] == "Y":
        print("Last Exon Error!!")
        assert False

    
    complement,scaffold = upstream_exon_line[4], upstream_exon_line[1]

    if complement == "0":
        gene_start = min(int(upstream_exon_line[2]),int(upstream_exon_line[3]))
        gene_end = max(int(downstream_exon_line[2]),int(downstream_exon_line[3]))
    
    if complement == "1":
        gene_start = min(int(downstream_exon_line[2]),int(downstream_exon_line[3]))
        gene_end = max(int(upstream_exon_line[2]),int(upstream_exon_line[3]))

    
    gene_sequence = get_gene_sequence(genome_location, species, genome_file, scaffold,gene_start,gene_end,complement,annotated_genome_location,annotated_species_name)
    print("BB")
    return(gene_sequence,scaffold,gene_start,gene_end )
    



def get_annotated_genome_name(annotated_genome_location, species):
    list_of_annotated_genomes = os.listdir(annotated_genome_location)
  
    # print(annotated_species)
    for annotated_species in list_of_annotated_genomes:
        if annotated_species.endswith(species):
            return(annotated_species)
             
    if annotated_species_name == '':
        print(f"Error with annotated species name")
        assert False
        
        


    

for species in species_list:


    annotated_genome_location = blast_output_location
    

    genome_file = get_genome_file(genome_location,species)
    annotated_species_name = get_annotated_genome_name(annotated_genome_location, species)

    gene_sequence,scaffold, gene_start, gene_end = get_sequence(blast_output_location,genome_location, species, genome_file,
                                                               annotated_genome_location,annotated_species_name)

    #print(gene_sequence,scaffold, gene_start, gene_end)
    with io.open(f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon/gene_sequence_all.fa",'w') as out_file:
        output = f">{species}_{scaffold}_{gene_start}_{gene_end}\n{gene_sequence}"
        out_file.write(output)

    local_genomic_fragment_location = f"{annotated_genome_location}/{annotated_species_name}/Period_gene_genomic_sequence_individual_exon/"
    cd_command = f'cd "{local_genomic_fragment_location}"\nmakeblastdb -in gene_sequence_all.fa -dbtype nucl\n'
        # # os.system(f'{cd_command}')
    subprocess.run(f'{cd_command}', shell = True, stderr = subprocess.DEVNULL)
        # # print(mkdb_command)
        # # os.system(f'{mkdb_command}')
    # print(gene_sequence)