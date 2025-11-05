import os

family_group_list = ["3.Satyrine","4.Pierinae","5.Coliadinae","6.Heliconiinae_Danainae_Nymphalinae","7.Papilionidae","8.Lycaenidae"]
location = "/mnt/griffin/saubar/Working_folder_circadian_rhythm/11.CWO_Exon_Analysis"

for family_group in family_group_list:
    list_of_species = os.listdir(f"{location}/{family_group}/1.Blast_result")
    if "desktop.ini" in list_of_species:
        list_of_species.remove("desktop.ini")
        
    for species in list_of_species:
        with open (f"{location}/{family_group}/1.Blast_result/{species}/final_coordinates.csv", 'r') as coordinate_file:
            coordinate_list = coordinate_file.readlines()
        
        scaffold = coordinate_list[1].split(",")[1]
        
        for i in range(2,len(coordinate_list)):
            coordinate_list[i]= coordinate_list[i].replace(",,", f",{scaffold},")
        
        print("".join(coordinate_list))
        output = ("".join(coordinate_list))
        with open (f"{location}/{family_group}/1.Blast_result/{species}/final_coordinates.csv", 'w') as coordinate_file:
            coordinate_file.write(output)
        
