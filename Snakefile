configfile: "config/config.yaml"
chroms=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]

import subprocess
import os
import pandas as pd

dbNSFP_file=pd.read_csv("config/dbnsfp_files.txt", names=["uniprot_ids"])
PDB_file=pd.read_csv("config/pdb_ids.txt", names=["PDB_ID"])
PDB_file[["prefix","uniprot_ids","model","postfix"]]=PDB_file["PDB_ID"].str.split("-", expand=True,)
PDB_dbNSFP=PDB_file[PDB_file.uniprot_ids.isin(dbNSFP_file.uniprot_ids)]
relevant_uniprot_ids=PDB_dbNSFP["uniprot_ids"].tolist()
relevant_uniprot_ids=list(set(relevant_uniprot_ids)) # remove duplicates
relevant_alphafold_models=PDB_dbNSFP["PDB_ID"].tolist()

grid_search_table=pd.read_csv(filepath_or_buffer="resources/grid_search.tsv", sep="\t").astype(str)
#grid_search_table=grid_search_table.iloc[[0,1]] # testing


rule all:
	input:
		#expand("data/pdb_features/{pdb}/b_factors.csv", pdb=relevant_alphafold_models),
		#expand("data/network/{pdb_name}_combined_w_network_and_dbnsfp.csv.gz", pdb_name=relevant_alphafold_models),
		#expand("data/split_dbNSFP/chr{chr}_ok",chr=chroms),
		"data/train_testset1/gnomad_extracted.csv.gz",
		#"data/merge_all/all_possible_values_concat.csv.gz",
		#expand("data/prediction/{prefix}_results.tsv", prefix=grid_search_table["prefix"].to_list()),
		#"data/joined_grid/joined_grid.tsv",
		"data/prediction/final_written_full_model.RData",
		"data/validation_set/validation_set_w_AlphScore.csv.gz",
		"data/analyse_score/spearman_plot.pdf",
		"data/plot_k/barplot_preprocessed.pdf",
		"data/combine_scores/aucs.pdf"


rule download_dbNSFP_AlphaFold_files:
	output:
		"data/dbNSFP_raw/download_ok",
		expand(config["pdb_dir"]+"{pdb_name}.pdb.gz", pdb_name=relevant_alphafold_models)
	params:
		dbNSFP="data/dbNSFP_raw",
		alphafold="data/pdb_files",
		partition=config["long_partition"]
	resources: time_job=4800, mem_mb=8000
	shell:
		"""
		homedir=$(pwd)
		mkdir -p "{params.dbNSFP}"
		mkdir -p "{params.alphafold}"
		
		cd "{params.alphafold}"
		rm -f UP000005640_9606_HUMAN_v2.tar
		wget {config[alphafold_download_adress]}
		tar -xf UP000005640_9606_HUMAN_v2.tar
		rm UP000005640_9606_HUMAN_v2.tar
		
		cd $homedir
		
		cd {params.dbNSFP}
		rm -f dbNSFP4.2a.zip
		wget {config[dbNSFP_download_adress]}
		unzip dbNSFP4.2a.zip
		rm dbNSFP4.2a.zip
		
		cd $homedir
		touch {output}
		"""


rule split_dbNSFP:
	input:
		"data/dbNSFP_raw/download_ok",
	output:
		"data/split_dbNSFP/chr{chr}_ok"
	params:
		tmp="data/split_dbNSFP/tmp/chr{chr}/",
		db_file="data/dbNSFP_raw/dbNSFP4.2a_variant.chr{chr}.gz",
		outdir="data/split_dbNSFP/by_uniprotID/",
		partition=config["short_partition"]
	resources: time_job=480, mem_mb=5000
	shell:
		"""
		homedir=$(pwd)
		mkdir -p "{params.tmp}"
		mkdir -p "{params.outdir}"
		cd "{params.tmp}"
		mkdir -p header
		
		if (( $(zcat $homedir"/{params.db_file}" | head -n1 > header/header.txt) )) # pipe / head combination produces non-0-exit
		then
		echo "not ok"
		fi
		
		zcat $homedir"/{params.db_file}" | tail -n+2 | split -l 20000 - {wildcards.chr}
		
		for file in $(find * -maxdepth 0 -type f)
		do
		python $homedir"/scripts/split_dbNSFP_unnest.py" $homedir"/config/pdb_ids.txt" $(pwd)"/header/header.txt" $file $homedir"/"{params.outdir}
		done
		
		cd $homedir
		touch "{output}"
		"""

rule get_FEATURE_DSSP:
	input:
		pdb=config["pdb_dir"]+"{pdb_name}.pdb.gz"
	output:
		"data/pdb_features/{pdb_name}/{pdb_name}_core.txt"
	params:
		partition=config["short_partition"],
	resources: time_job=20
	shell:
		"""
		feature="$(pwd)/"{config[feature_folder]}"bin/featurize"
		export FEATURE_DIR="$(pwd)/"{config[feature_folder]}"data/"
		work_dir="data/pdb_features/{wildcards.pdb_name}/"
		temp_pdb=$work_dir"1A1A.pdb"
		
		export DSSP_DIR=$work_dir
		export PDB_DIR=$work_dir
		
		mkdir -p $work_dir
		gunzip -c "{input}" > $temp_pdb
		
		export LD_LIBRARY_PATH={config[dssp_dir]}
		
		#run dssp
		cat $temp_pdb | {config[dssp_bin]} -i /dev/stdin > $work_dir"1A1A.dssp"
		echo 1A1A > $work_dir"pdb_id"
		echo  $work_dir"1A1A.dssp"
		
		$feature $temp_pdb -n1 -w0.2 | grep ":A@CA$" > $work_dir"{wildcards.pdb_name}_core.txt"
		$feature $temp_pdb -n1 -w3 | grep ":A@CA$" > $work_dir"{wildcards.pdb_name}_3A.txt"
		$feature $temp_pdb -n1 -w6 | grep ":A@CA$" > $work_dir"{wildcards.pdb_name}_6A.txt"
		$feature $temp_pdb -n1 -w9 | grep ":A@CA$" > $work_dir"{wildcards.pdb_name}_9A.txt"
		$feature $temp_pdb -n1 -w12 | grep ":A@CA$" > $work_dir"{wildcards.pdb_name}_12A.txt"
		
		rm $temp_pdb*
		"""

rule get_feature_protinter_interactions:
	input:
		pdb=config["pdb_dir"]+"{pdb_name}.pdb.gz"
	output:
		"data/pdb_features/{pdb_name}/result_hbond_side_side.csv"
	params:
		partition=config["short_partition"]
	resources: time_job=60
	shell:
		"""
		protinter="$(pwd)/"{config[protinter_bin]}
		work_dir="data/pdb_features/{wildcards.pdb_name}/"
		mkdir -p $work_dir
		
		gunzip -c "{input}" > $work_dir"/pdbfile.pdb"
		cd $work_dir
		
		python3 $protinter pdbfile.pdb -hydrophobic -csv
		python3 $protinter pdbfile.pdb -disulphide -csv
		python3 $protinter pdbfile.pdb -ionic -csv
		python3 $protinter pdbfile.pdb -aroaro -csv
		python3 $protinter pdbfile.pdb -arosul -csv
		python3 $protinter pdbfile.pdb -catpi -csv
		#python3 $protinter pdbfile.pdb -hb1 -csv
		python3 $protinter pdbfile.pdb -hb2 -csv
		python3 $protinter pdbfile.pdb -hb3 -csv
		python3 $protinter pdbfile.pdb -within_radius -csv
		
		sed -i "s/ //g" result_*
		rm pdbfile.pdb
		"""

rule get_feature_HSE:
	input:
		pdb=config["pdb_dir"]+"{pdb_name}.pdb.gz"
	output:
		"data/pdb_features/{pdb_name}/HSEs.csv"
	params:
		pdb_path=config["pdb_dir"]+"{pdb_name}.pdb",
		partition=config["short_partition"]
	resources: time_job=20
	run:
		import os
		from Bio.PDB.PDBParser import PDBParser
		import Bio.PDB as bpd
		import pandas as pd
		pdb_path=params[0]
		os.system('gunzip -c "' + pdb_path + '.gz" > "' + pdb_path + 'X"' )
		
		parser = PDBParser()
		structure = parser.get_structure("test", pdb_path + 'X')
		hse = bpd.HSExposure
		exp_ca = hse.HSExposureCB(structure)
		os.system('rm "' + pdb_path + 'X"' )
		resi_count=len(exp_ca.keys())
		print("ja")
		list_HSE=[ [exp_ca.keys()[i][1][1], exp_ca.property_list[i][1][0], exp_ca.property_list[i][1][1] ] for i in range(0, resi_count)]
		df_HSE=pd.DataFrame(list_HSE)
		df_HSE.columns =['RESI','HSE1','HSE2']
		df_HSE.index=df_HSE.RESI
		df_HSE.to_csv(output[0], index=False)


rule get_feature_pLDDT:
	input:
		pdb=config["pdb_dir"]+"{pdb_name}.pdb.gz"
	output:
		"data/pdb_features/{pdb_name}/b_factors.csv"
	params:
		partition=config["short_partition"]
	resources: time_job=20
	run:
		import pandas as pd 
		from biopandas.pdb import PandasPdb
		import numpy as np
		ppdb_load = PandasPdb()
		ppdb_load.read_pdb(input[0])
		pd_df=ppdb_load.df["ATOM"][["residue_number","atom_name","b_factor"]].to_numpy()
		np.savetxt(output[0], pd_df[pd_df[:,1]=="CA"], fmt='%s')


rule combine_features_pdb_level:
	input:
		"data/pdb_features/{pdb_name}/result_hbond_side_side.csv",
		"data/pdb_features/{pdb_name}/b_factors.csv",
		"data/pdb_features/{pdb_name}/{pdb_name}_core.txt",
		"data/pdb_features/{pdb_name}/HSEs.csv"
	output:
		"data/combine1_pdb_data/{pdb_name}_combined_features.csv.gz"
	params:
		pdb_directory="data/pdb_features/{pdb_name}/",
		header_feature="resources/header_featurize.txt",
		partition=config["short_partition"]
	resources: time_job=480, mem_mb=8000		
	script:
		"scripts/combine_features_pdb_level.py"


rule add_graph_based_analysis_pdb_level:
	input:
		comb1="data/combine1_pdb_data/{pdb_name}_combined_features.csv.gz",
		pdb=config["pdb_dir"]+"{pdb_name}.pdb.gz",
		dbnsfp=expand("data/split_dbNSFP/chr{chr}_ok",chr=chroms),
	output:
		"data/network/{pdb_name}_combined_w_network_and_dbnsfp.csv.gz"
	params:
		"data/split_dbNSFP/by_uniprotID/",
		pdb_path=config["pdb_dir"]+"{pdb_name}.pdb",
		partition=config["short_partition"]
	resources: time_job=240, mem_mb=15000
	script:
		"scripts/graph_based_analysis_pdb_level.py"
		

rule combine_pdb_level_files_to_protein_level:
	input:
		lambda wildcards : ["data/network/" + el.rstrip("\n") + "_combined_w_network_and_dbnsfp.csv.gz" for el in relevant_alphafold_models if el.split("-")[1] in wildcards.uniprot_id] # AF-A0A087WUL8-F3-model_v1_combined_features
	output:
		"data/combine2_protein/{uniprot_id}_w_AFfeatures.csv.gz"
	params:
		"data/split_dbNSFP/by_uniprotID/",
		partition=config["short_partition"]
	resources: time_job=240, mem_mb=lambda wildcards : 3000+1000*len([el.rstrip("\n") for el in relevant_alphafold_models if el.split("-")[1] in wildcards.uniprot_id])
	run:
		import pandas as pd
		import os
		li =[]
		for filename in input:
			try:
				df = pd.read_csv(filename, low_memory=False)
				li.append(df)
			except:
				print(filename + " error")
		try:
			frame = pd.concat(li, axis=0, ignore_index=True)
		except:
			print("empty df")
			os.system('touch ' + output[0])
        		os._exit(0)
		dbnsfp_cols=pd.read_csv(params[0] + wildcards[0] + ".csv.gz", na_values=".", low_memory=False).columns.to_list()
		numeric_cols=frame.select_dtypes(include=['number']).columns.to_list()
		numeric_cols_not_in_dbnsfp=[col for col in numeric_cols if col not in dbnsfp_cols]
		median_vals=frame.groupby(['pos(1-based)','ref','alt'])[numeric_cols_not_in_dbnsfp].median()
		constant_vals=frame.drop_duplicates(subset=['pos(1-based)','ref','alt'])
		constant_vals=constant_vals.drop(numeric_cols_not_in_dbnsfp, axis=1)
		new_vals=pd.merge(constant_vals, median_vals,  how='left', left_on=['pos(1-based)','ref','alt'], right_on = ['pos(1-based)','ref','alt'])
		new_vals["protein_mean_CADD"]=new_vals["protein_mean_CADD"].median() 
		new_vals["protein_mean_b_factor"]=new_vals["protein_mean_b_factor"].median()
		new_vals["protein_length"]=new_vals["RESI"].max() 
		new_vals=new_vals.sort_values(by=["RESI"], axis=0)
		new_vals.to_csv(output[0])


rule extract_clinvar_gnomad_sets:
	input:
		expand("data/combine2_protein/{uniprot_id}_w_AFfeatures.csv.gz", uniprot_id=relevant_uniprot_ids)
	output:
		csv_gz_out="data/train_testset1/gnomad_extracted.csv.gz"
	params:
		partition=config["short_partition"],
		loaction_csvgzs="data/combine2_protein/"
	resources: time_job=480, mem_mb=64000
	script:
		"scripts/extract_clinvar_gnomad_sets.R"

rule preprocess_clinvar_gnomad_set:
	input:
		"data/train_testset1/gnomad_extracted.csv.gz"
	output:
		preprocessed="data/train_testset1/gnomad_extracted_prepro.csv.gz",
		recalibrated="data/train_testset1/gnomad_extracted_prepro_rec.csv.gz"
	resources: cpus=1, mem_mb=18000, time_job=480
	params:
		partition=config["short_partition"]
	shell:
		"""
		Rscript scripts/preprocess.R -i {input} -o {output.preprocessed}
		Rscript scripts/recalibrate_variants.R -i {output.preprocessed} -o {output.recalibrated}
		"""

rule get_properties_clinvar_gnomad_set:
	input:
		preprocessed="data/train_testset1/gnomad_extracted_prepro.csv.gz",
		recalibrated="data/train_testset1/gnomad_extracted_prepro_rec.csv.gz",
		excel="resources/available_colnames_W_surr.xlsx"
	output:
		barplot_pre="data/plot_k/barplot_preprocessed.pdf",
		barplot_rec="data/plot_k/barplot_recalibrated.pdf",
		AA_comb="data/plot_k/AA_ex_compared.pdf",
		char_var="data/plot_k/characteristics_variants.tsv"
	resources: cpus=1, mem_mb=18000, time_job=480
	params:
		partition=config["short_partition"],
		out_folder="data/plot_k/"
	shell:
		"""
		Rscript scripts/Fig_recal.R \
		--input_preprocessed {input.preprocessed} \
		--input_recalibrated {input.recalibrated} \
		--out_folder {params.out_folder} 

		Rscript scripts/AA_freq.R \
		--input_recalibrated {input.recalibrated} \
		--excel_location {input.excel} \
		--out_folder {params.out_folder} 

		Rscript scripts/char_variants.R \
		--csv_location {input.recalibrated} \
		--out_folder {params.out_folder} 
		"""

rule training_do_grid_search:
	input:
		"resources/available_colnames_W_surr.xlsx",
		"data/train_testset1/gnomad_extracted_prepro_rec.csv.gz"
	output:
		"data/prediction/{prefix}_results.tsv"
	resources: cpus=6, mem_mb=24000, time_job=480
	params:
		partition=config["short_partition"]
	run:
		import os
		parameters=grid_search_table[grid_search_table["prefix"]==wildcards[0]]
		os.system('Rscript scripts/prediction.R --prefix ' + 
			wildcards[0] + " --csv_location " + input[1] + 
			" " + " ".join(parameters.iloc[0,1:].to_list()))

rule join_grid_search_files:
	input:
		expand("data/prediction/{prefix}_results.tsv", prefix=grid_search_table["prefix"].to_list()),
	output:
		"data/joined_grid/joined_grid.tsv"
	resources: cpus=1, mem_mb=18000, time_job=480
	params:
		partition=config["short_partition"]
	shell:
		"""
		mkdir -p "data/joined_grid"
		cat {input} | tail -n3 | grep "auc_Alph_train_gnomAD" > {output}
		cat {input} | grep -v "auc_Alph_train_gnomAD" >> {output}
		"""

rule fit_models_w_final_settings_from_grid_search:
	input:
		excel="resources/available_colnames_W_surr.xlsx",
		csv="data/train_testset1/gnomad_extracted_prepro_rec.csv.gz"
	output:
		"data/prediction/final_written_full_model.RData",
		"data/prediction/final_toAS_properties.RData",
		"data/prediction/final_colnames_to_use.RData",
		"data/prediction/pre_final_model_variants.csv.gz"
	resources: cpus=8, mem_mb=30000, time_job=480
	params:
		partition=config["short_partition"]
	shell:
		"""
		Rscript scripts/prediction.R \
		--prefix pre_final_model_k_fold \
		--csv_location {input.csv} \
		--excel_location {input.excel} \
		--method_pred extratree \
		--num_trees_param 2000 \
		--max_depth_param 6 \
		--min_node_param 10 \
		--cor_param 0.999999 \
		--b_factor_param 75 \
		--eta_param 0 \
		--subsample_param 0 \
		--min_child_weight_param 0 \
		--gamma_param 0 \
		--k_fold_cross_val TRUE

		Rscript scripts/prediction.R \
		--prefix pre_final_model \
		--csv_location {input.csv} \
		--excel_location {input.excel} \
		--method_pred extratree \
		--num_trees_param 2000 \
		--max_depth_param 6 \
		--min_node_param 10 \
		--cor_param 0.999999 \
		--b_factor_param 75 \
		--eta_param 0 \
		--subsample_param 0 \
		--min_child_weight_param 0 \
		--gamma_param 0 \
		--write_model TRUE \
		--write_dataset TRUE \
		--importance impurity
		
		Rscript scripts/prediction.R \
		--prefix pre_final_model \
		--csv_location {input.csv} \
		--excel_location {input.excel} \
		--method_pred extratree \
		--num_trees_param 2000 \
		--max_depth_param 6 \
		--min_node_param 10 \
		--cor_param 0.999999 \
		--b_factor_param 75 \
		--eta_param 0 \
		--subsample_param 0 \
		--min_child_weight_param 0 \
		--gamma_param 0 \
		--write_model FALSE \
		--write_dataset FALSE \
		--importance permutation

		Rscript scripts/prediction.R \
		--prefix final \
		--csv_location {input.csv} \
		--excel_location {input.excel} \
		--method_pred extratree \
		--num_trees_param 2000 \
		--max_depth_param 6 \
		--min_node_param 10 \
		--cor_param 0.999999 \
		--b_factor_param 75 \
		--eta_param 0 \
		--subsample_param 0 \
		--min_child_weight_param 0 \
		--gamma_param 0 \
		--full_model TRUE \
		--write_model TRUE
		"""

rule predict_Alphscore_protein_level:
	input:
		csv="data/combine2_protein/{uniprot_id}_w_AFfeatures.csv.gz",
		model="data/prediction/final_written_full_model.RData",
		to_AS="data/prediction/final_toAS_properties.RData",
		colnames="data/prediction/final_colnames_to_use.RData",
	output:
		"data/predicted_prots/{uniprot_id}_w_AlphScore_red_{reduced}.csv.gz"
	resources: cpus=1, time_job=60, mem_mb=lambda wildcards : 4000+1000*len([el.rstrip("\n") for el in relevant_alphafold_models if el.split("-")[1] in wildcards.uniprot_id])
	params:
		partition=config["short_partition"]
	shell:
		"""
		Rscript scripts/predict_w_ranger_model.R \
		--csv_location {input.csv} \
		--model_location {input.model} \
		--output_file {output} \
		--use_cols_file {input.colnames} \
		--toAS_properties {input.to_AS} \
		--reduced {wildcards.reduced}
		"""

rule create_validation_set_DMS:
	input:
		"resources/",
		expand("data/predicted_prots/{uniprot_id}_w_AlphScore_red_FALSE.csv.gz",
		 uniprot_id=["P01112","P04637","P05067",
		"P07550","P0DP23","P38398","P43246","P51580","P60484","P63165","P63279","Q9BQB6","Q9H3S4"])
	output:
		compiled="data/validation_set/validation_set_w_AlphScore.csv.gz",
	resources: cpus=1, mem_mb=18000, time_job=480
	params:
		partition=config["short_partition"]
	shell:
		"""
		Rscript scripts/compile_scores.R
		"""

rule analyse_performance_on_DMS_data:
	input:
		compiled="data/validation_set/validation_set_w_AlphScore.csv.gz",
		test_dataset="data/prediction/pre_final_model_variants.csv.gz"
	output:
		one_plot="data/analyse_score/spearman_plot.pdf",
	resources: cpus=1, mem_mb=18000, time_job=480
	params:
		partition=config["short_partition"],
		out_folder="data/analyse_score/"
	shell:
		"""
		Rscript scripts/analyse_score.R \
		--variants {input.test_dataset} \
		--validation_set {input.compiled} \
		--out_folder {params.out_folder}
		"""

rule combine_alphafold_w_existing_scores:
	input:
		variant_dataset="data/prediction/pre_final_model_variants.csv.gz"
	output:
		one_plot="data/combine_scores/aucs.pdf",
	resources: cpus=1, mem_mb=18000, time_job=480
	params:
		partition=config["short_partition"],
		out_folder="data/combine_scores/"
	shell:
		"""
		Rscript scripts/combine_scores.R \
		--variants {input.variant_dataset} \
		--out_folder {params.out_folder}
		"""

rule combine_proteins_level_w_Alphscore_to_one_file:
	input:
		expand("data/predicted_prots/{uniprot_id}_w_AlphScore_red_TRUE.csv.gz", uniprot_id=relevant_uniprot_ids)
	output:
		"data/merge_all/all_possible_values_concat.csv.gz"
	params:
		partition=config["long_partition"],
		in_folder="data/predicted_prots/",
		out_folder="data/merge_all/"
	resources: time_job=4800, mem_mb=8000
	shell:
		"scripts/combine_all_proteins_to_one_file.sh {params.in_folder} {params.out_folder}"
