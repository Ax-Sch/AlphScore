import sys
import pandas as pd
import numpy as np
from operator import itemgetter


### for debuging
class testclass:
	def __init__(self):
		self.wildcards=["AF-A0A087WUL8-F11-model_v1"]
		self.output=["data/AF-A0A087WUL8-F11-model_v1__combined_features.csv.gz"]
		self.params=["data/pdb_features/"]

#snakemake=testclass()
offset=0
offset=( int(snakemake.wildcards[0].split("-")[2][1:] ) - 1) * 200

# generate a list with the number of interactions
possible_interactions=["aroaro", "arosul", "cationpi", "disulphide", "hbond_main_side", "hbond_side_side", "hydrophobic", "ionic", "within_radius"]
interaction_matrix=np.zeros((40000,len(possible_interactions)),np.int) # 40000 - longer than Titin
i=0
edges=[]
for interaction in possible_interactions:
	interactions_pd=pd.read_csv(snakemake.params[0] + "/result_" + interaction + ".csv", low_memory=False)
	interactions_np=interactions_pd.to_numpy()
	interactions_np[:,1] = interactions_np[:,1] + offset
	interactions_np[:,3] = interactions_np[:,3] + offset
	interactions_joined=np.append(interactions_np[:,1], interactions_np[:,3])
	interactions_np[:,0]=np.array([s.strip() for s in interactions_np[:,0]])
	interactions_np[:,2]=np.array([s.strip() for s in interactions_np[:,2]])
	for ele in interactions_joined:
		interaction_matrix[ele,i]=interaction_matrix[ele,i]+1
	i=i+1

		
# Read features from FEATURIZE Framework
header=pd.read_csv(snakemake.params[1],skipinitialspace=True).columns.to_list()

to_remove=["RAUTE","coord1", "coord2", "coord3", "RAUTE2", "SECONDARY_STRUCTURE2_IS_UNKNOWN", "SECONDARY_STRUCTURE1_IS_UNKNOWN", "RESIDUE_CLASS2_IS_UNKNOWN","RESIDUE_CLASS1_IS_UNKNOWN", "RESIDUE_NAME_IS_OTHER", "ATOM_TYPE_IS_OTHER", "ELEMENT_IS_OTHER", "Atom_name"]

for suff in ["_core", "_3A", "_6A", "_9A", "_12A"]:
	header_surround=[ el + suff  if el not in ["RESN_RESI_AT"] + to_remove else el for el in header] 
	features=pd.read_csv(snakemake.params[0] + "/" + snakemake.wildcards[0] + suff + ".txt", delimiter="\t", names=header_surround,skipinitialspace=True , low_memory=False)
	features["RESN_RESI"]=features.RESN_RESI_AT.str.split(":",expand=True).iloc[:,0]
	features["RESN"]=features.RESN_RESI.str[:3]
	features["RESI"]=features.RESN_RESI.str[3:]
	features=features.astype({"RESI": int})
	features.RESI=features.RESI + offset
	features["RESN_RESI"]=features["RESN"]+features["RESI"].astype(str)
	features.index=features.RESI
	features=features.drop(to_remove, axis=1)
	print(suff)
	if suff!="_core":
		features_joined=features_joined.join(features.drop(['RESN_RESI_AT', 'RESN_RESI', 'RESI','RESN'], axis=1))
	else:
		features_joined=features

# join interactions and FEATURE based values
interactions_for_join=pd.DataFrame(interaction_matrix, columns=possible_interactions) # offset already inserted
interactions_for_join=interactions_for_join.astype(int)
interactions_for_join.index.name="RESI" 
features_joined2=features_joined.join(interactions_for_join)

# join the bfactors
b_factors=pd.read_csv(snakemake.params[0] + "/b_factors.csv", delimiter=" ", names=["residue_number","atom_name","b_factor"], skipinitialspace=True, low_memory=False)
b_factors["residue_number"] = pd.to_numeric(b_factors["residue_number"]) + offset
b_factors.index=b_factors["residue_number"] ########### +
features_joined3=features_joined2.join(b_factors)

# join the HSEs
HSE=pd.read_csv(snakemake.params[0] + "/HSEs.csv", delimiter=",", skipinitialspace=True, low_memory=False)
print(HSE.columns)
HSE["RESI"] = pd.to_numeric(HSE["RESI"]) + offset
HSE.index=HSE["RESI"]
HSE=HSE.drop("RESI", axis=1) ########### +
features_joined4=features_joined3.join(HSE)


features_joined4[["pdb_file"]]=snakemake.wildcards[0]
features_joined4.to_csv(snakemake.output[0])
