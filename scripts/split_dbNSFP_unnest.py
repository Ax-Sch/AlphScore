# Unnest the columns Uniprot_acc_split and HGVSp_VEP_split in dbNSFP which contain multiple values (separated by ";")
# argument 1: file containing filenames of Alphafold2 pdbs
# argument 2: file containing the column-names of the dbNSFP file used
# argument 3: dbNSFP file to read and process
# argument 4: folder to save files to

import glob
import pandas as pd
import subprocess
import sys

print(sys.argv)

# get the relevant Uniprot-IDs from the file given by argument 1
my_file = open(sys.argv[1], "r")
pdb_ids = my_file.readlines()
pdb_ids2=[]
for element in pdb_ids:
    pdb_ids2.append(element.strip())
Uniprot_ids_in_alphafold=[line.split("-")[1] for line in pdb_ids2]
Uniprot_ids_in_alphafold=set(Uniprot_ids_in_alphafold)

# read argument 2 / 3, the header and the file containing the chunk of dbNSFP 
header_file = pd.read_csv(sys.argv[2], sep="\t")
dbNSFP_split_part=pd.read_csv(sys.argv[3], sep="\t", names=header_file.columns.to_list())

# convert the fields in Uniprot_acc_split / HGVSp_VEP_split to lists
dbNSFP_split_part["Uniprot_acc_split"]=dbNSFP_split_part.Uniprot_acc.str.split(";")
dbNSFP_split_part["HGVSp_VEP_split"]=dbNSFP_split_part.HGVSp_VEP.str.split(";")

# unnest dbNSFP_split_part.HGVSp_VEP_split e.g. from [["a1","a2","a3"],["b1","b2"], ["c1"]] to ["a1","a2","a3","b1","b2","c1"]
hgvs_unnested=[item for sublist in dbNSFP_split_part.HGVSp_VEP_split for item in sublist]
Uniprot_unnested=[item for sublist in dbNSFP_split_part.Uniprot_acc_split for item in sublist]

# create an index "index_list" that e.g. looks like this [1 1 1 2 2 3] if dbNSFP_split_part.Uniprot_acc_split contained [["a1","a2","a3"],["b1","b2"],["c1"]] 
lengths_sublists=([len(sublist) for sublist in dbNSFP_split_part.Uniprot_acc_split ])
i=0
k=[]
for i in range(0,len(dbNSFP_split_part)):
	k.append([i] * lengths_sublists[i])
	i=i+1
index_list = [item for sublist in k for item in sublist]

# if a row in dbNSFP_split_part contains n entries in HGVSp_VEP, then this row will be copied n times. Then each of the n rows can be attributed to a single value from HGVSp_VEP.
dbNSFP_split_part=dbNSFP_split_part.iloc[index_list,:]
dbNSFP_split_part.Uniprot_acc_split=Uniprot_unnested
dbNSFP_split_part.HGVSp_VEP_split=hgvs_unnested

# filter for relevant Uniprot IDs
dbNSFP_split_part2=dbNSFP_split_part[dbNSFP_split_part.Uniprot_acc_split.isin(Uniprot_ids_in_alphafold)]

# save to files. A split could have occured in the middle of a gene - therefore a file for a given gene / uniprot-ID could already be present, which will cause an exception and is handled by the except part.
for unip_id in set(dbNSFP_split_part2.Uniprot_acc_split.to_list()):
	file_base=sys.argv[4] + unip_id
	try:
		dbNSFP_split_part2[dbNSFP_split_part2.Uniprot_acc_split==unip_id].to_csv(file_base + ".csv.gz", 
		   index=False, 
		   compression="gzip", mode="x")
	except:
		dbNSFP_split_part2[dbNSFP_split_part2.Uniprot_acc_split==unip_id].to_csv(file_base + "_2.csv.gz", 
		   index=False, compression="gzip", mode="x", header=False)
		bash_command = "zcat "+ file_base + ".csv.gz " +  file_base + "_2.csv.gz | bgzip > " + file_base + "temp" # combine to temp file
		print(bash_command)
		subprocess.run(bash_command, shell=True)
		bash_command = "mv " + file_base + "temp" + " " + file_base + ".csv.gz"
		print(bash_command)
		subprocess.run(bash_command, shell=True)
		bash_command= "rm -f " + file_base + "_2.csv.gz"
		print(bash_command)
		subprocess.run(bash_command, shell=True)
