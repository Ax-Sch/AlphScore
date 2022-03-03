import sys
import pandas as pd
import numpy as np
from operator import itemgetter
import networkx as nx
import biographs as bg
import os

def exit_procedure(text_to_log):
	os.system('touch ' + snakemake.output[0])
	outF = open(snakemake.output[0]+".log", "w")
	outF.write(text_to_log)
	outF.close()
	print(text_to_log)
	os._exit(0)

#results/{pdb_name}_combined_features.csv.gz
class testclass:
	def __init__(self):
		self.wildcards=["AF-A0A087WUL8-F11-model_v1"]
		self.output=["results/AF-A0A087WUL8-F11-model_v1_features.tsvneighbour_vals.csv"]
		self.input=["results/AF-A0A087WUL8-F11-model_v1__combined_features.csv.gz"]
		self.params=["/media/axel/Dateien/Arbeit_Gen/alphafold2/dbnsfp_results/"]

#snakemake=testclass()

offset=0
offset=( int(snakemake.wildcards[0].split("-")[2][1:] ) - 1) * 200

#print(snakemake.wildcards, snakemake.output, snakemake.params)


pdb_path=snakemake.params[1]
os.system('gunzip -c "' + pdb_path + '".gz > ' +pdb_path)
molecule = bg.Pmolecule(pdb_path)
os.system('rm -f "' + pdb_path+ '"')

# biopython molecule structural model
mol_model = molecule.model
# networkx graph, by default 5 angstrom
network = molecule.network()
network = molecule.network(cutoff=4, weight=True)

start_resi=int(list(network.nodes)[0][1:])
network_relabeled=nx.convert_node_labels_to_integers(network, start_resi + offset)



#### read combined features

features_joined3=pd.read_csv(snakemake.input[0] ,low_memory=False)

#### read dbNSFP

dbnsfp=pd.read_csv(snakemake.params[0] + snakemake.wildcards[0].split("-")[1] + ".csv.gz", na_values=".", low_memory=False)


#### aggregate features of nodes
cons_scores=['CADD_raw','CADD_phred','BayesDel_noAF_score','REVEL_score', 'GERP++_NR', 'GERP++_RS', 'GERP++_RS_rankscore', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP30way_mammalian', 'phyloP30way_mammalian_rankscore', 'phyloP17way_primate', 'phyloP17way_primate_rankscore', 'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons30way_mammalian', 'phastCons30way_mammalian_rankscore', 'phastCons17way_primate', 'phastCons17way_primate_rankscore'] # removed:  'Polyphen2_HDIV_score', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'SIFT_score', 'SIFT_converted_rankscore', 'SIFT_pred', 'SiPhy_29way_pi'

def convert_int(num_str):
	try:
		num_int=int(num_str)
		return(num_int)
	except:
		return(0)

dbnsfp.HGVSp_VEP_split=dbnsfp.HGVSp_VEP_split.replace("p.Met1?","p.Met1???")
num_string=dbnsfp.HGVSp_VEP_split.str[5:-3]
dbnsfp["RESI"]=[convert_int(item) for item in num_string.to_list()]
dbnsfp["RESN"]=dbnsfp.HGVSp_VEP_split.str[2:5]
dbnsfp.RESN=dbnsfp.RESN.str.upper()
dbnsfp["RESN_RESI"]=dbnsfp["RESN"]  + dbnsfp["RESI"].astype(str)
dbnsfp=dbnsfp[dbnsfp.RESI!=0]
dbnsfp["to_AS"]=dbnsfp.HGVSp_VEP_split.str[-3:]
dbnsfp=dbnsfp[dbnsfp.to_AS!="Ter"]
dbnsfp_inf=dbnsfp.groupby('RESI')[cons_scores].mean()


protein_mean_CADD=dbnsfp_inf.mean()["CADD_raw"]
protein_mean_b_factor=features_joined3["b_factor"].mean()

dbnsfp_resn_resi=set(dbnsfp.RESN_RESI)

mismatch_count=0
nodes=[]

resis=list(set(features_joined3["RESI"].tolist()))

for i in resis:
	try:
		node_dict=features_joined3[features_joined3['RESI']==i].to_dict('records')[0]
		#if i in dbnsfp_inf.index:
		if not (node_dict["RESN_RESI"] in dbnsfp_resn_resi):
			print(node_dict["RESN_RESI"], "not present in dbnsfp!")
			mismatch_count=mismatch_count+1
			#import pdb; pdb.set_trace()
		node_dict={**node_dict, **dbnsfp_inf[dbnsfp_inf.index==i].to_dict('records')[0]}
		nodes.append([i, node_dict])
		#import pdb; pdb.set_trace()
	except:
		exit_procedure("1 Node construction failed")


if mismatch_count>i/2:
	exit_procedure("2 mismatched amino acids: " + str(mismatch_count))


network_relabeled.add_nodes_from(nodes)
print("Network okay!")

def create_pd_dict(diction, name_):
	return(pd.DataFrame.from_dict(diction, orient="index", columns=[name_]))

#score="HYDROPHOBICITY_3A"
physicochemical_scores=["CHARGE_3A","HYDROPHOBICITY_3A","MOBILITY_3A","SOLVENT_ACCESSIBILITY_core",'HSE1','HSE2']
dict_vals={}
list_vals=pd.DataFrame()
for score in physicochemical_scores + cons_scores:
	print(score + " ok.")
	for x in range(int(list(network_relabeled.nodes)[0]),int(list(network_relabeled.nodes)[-1])+1):
		sum_weight=0
		sum_value=0
		for neighbor in network_relabeled.neighbors(x):
			values=network_relabeled.nodes[neighbor][score]
			weight=network_relabeled[x][neighbor]["weight"]
			sum_weight=sum_weight+weight
			sum_value=sum_value+weight*values
			#print(values, weight, neighbor)
		#print(dict_vals)
		dict_vals={**dict_vals, **{x :  sum_value/sum_weight}}
	if list_vals.size==0:
		df_join=create_pd_dict(dict_vals, score+"_sur")
		list_vals=df_join.fillna(df_join.median())
	else:
		df_join=create_pd_dict(dict_vals, score+"_sur")
		list_vals=list_vals.join(df_join.fillna(df_join.median()))



deg=create_pd_dict(dict(nx.degree(network_relabeled)), "degree")
betw_cent=create_pd_dict(nx.betweenness_centrality(network_relabeled), "betweenness_centrality")
av_neigh=create_pd_dict(nx.average_neighbor_degree(network_relabeled), "av_neighbor_degree")
ev_cent=create_pd_dict(nx.eigenvector_centrality(network_relabeled, max_iter=100000, tol=0.0001),"eigenvec_centrality")

list_vals=list_vals.join([av_neigh, ev_cent, betw_cent, deg])
features_joined3.index = features_joined3.RESI

dbnsfp=dbnsfp.join(on="RESI", other=features_joined3, rsuffix="_features")
dbnsfp=dbnsfp.join(on="RESI", other=list_vals)
dbnsfp=dbnsfp[dbnsfp.residue_number.notnull()]

dbnsfp["protein_mean_CADD"]=protein_mean_CADD
dbnsfp["protein_mean_b_factor"]=protein_mean_b_factor

dbnsfp.to_csv(snakemake.output[0])

