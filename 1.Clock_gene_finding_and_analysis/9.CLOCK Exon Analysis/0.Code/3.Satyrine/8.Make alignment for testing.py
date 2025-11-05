from Bio import SeqIO
import os
import subprocess

file_location = os.getcwd()
family_group = file_location.split("/")[-1]

number_of_exons = 16

location = f"/mnt/h/My Drive/Circadian Rhythm Genes Project/7.Timeless Exon Analysis/{family_group}"

list_of_folders = os.listdir(location)
if "3.Test_alignment" in list_of_folders:
    subprocess.run(f'rm -r "{location}/3.Test_alignment"', shell = True, stderr = subprocess.DEVNULL)
os.mkdir(f"{location}/3.Test_alignment")

    
exons_location = f"{location}/2.Final_output"

list_of_species = os.listdir(exons_location)

if "desktop.ini" in list_of_species:
    list_of_species.remove("desktop.ini")


for exon_number in range(1,number_of_exons+1):

    exon_align = ''
    for species in list_of_species:
#         print(species, "Exon", exon_number)
        try:
            sequence = SeqIO.parse(f"{exons_location}/{species}/Exon_{exon_number}.fa", 'fasta')

            for records in sequence:

                left_overhang = int(records.id.split("_")[-3])

                sequence = records.seq[left_overhang:].translate()

                exon_align += f">{records.id}\n{sequence}\n\n"

        except:
            
            pass
    
    
    with open(f"{location}/3.Test_alignment/Exon_{exon_number}.fa",'w') as out_file:
        out_file.write(exon_align)


# Align using MAFFT
from Bio import SeqIO

error_exon_file = ''

list_of_files_to_run_mafft_on = os.listdir(f"{location}/3.Test_alignment/")
for file in list_of_files_to_run_mafft_on:
    print(file)
    if file.endswith(".fa"):
        command = f'"mafft" --localpair --clustalout --maxiterate 16 --reorder "{location}/3.Test_alignment/{file}" > "{location}/3.Test_alignment/alignment_clustal_{file}.txt"'
        subprocess.run(f'{command}', shell=True, stderr = subprocess.DEVNULL) 
    
        with open(f"{location}/3.Test_alignment/alignment_clustal_{file}.txt", 'r') as alignment_file:
            print("".join(alignment_file.readlines()))

        error_exon = input("continue?")
        if error_exon != "y":

            if "," not in error_exon:
                error_exon_list = [error_exon]
            else:
                error_exon_list = error_exon.split(",")

            for species in error_exon_list:
                sequence = SeqIO.parse(f"{location}/3.Test_alignment/{file}", "fasta")

                for records in sequence:
                    if species in records.id:
                        print(records.id)
                        species_name = input("Species_name?")
                        scaffold = input("Scaffold:")
                        exon = file.split("_")[1].split(".")[0]

                        error_exon_file += f"{species_name},{scaffold},000,000,0,Y,Error_Exon_{exon},00,00,00\n\n"


with open(f"{location}/3.Test_alignment/error_exons.txt", 'w') as out_file:
    out_file.write(error_exon_file)