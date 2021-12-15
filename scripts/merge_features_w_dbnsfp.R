library(data.table)
args <- commandArgs(trailingOnly = TRUE)

#testing:
#args[1]<-"/media/axel/Dateien/Arbeit_Gen/alphafold2/dbnsfp_data/Q8NGP4.csv.gz"
#args[2]<-"/media/axel/Dateien/Arbeit_Gen/alphafold2/extract_features_alphafold/data/AF-Q8NGP4-F1-model_v1_features.tsvneighbour_vals.csv"
#args[3]<-"/media/axel/Dateien/Arbeit_Gen/alphafold2/test.csv.gz"

dbnsfp<-fread(args[1])
dbnsfp[HGVSp_VEP_split=="p.Met1?", HGVSp_VEP_split:="p.Met1???"]
dbnsfp[,RESI:=as.numeric(substr(HGVSp_VEP_split,6,nchar(HGVSp_VEP_split)-3))]
dbnsfp[,from_AS:=substr(HGVSp_VEP_split,3,5)]
dbnsfp[,key_var:=paste(Uniprot_acc_split, RESI, from_AS, sep=":")]

features<-fread(args[2])
uniprot_id=strsplit(features$pdb_file, "-")[[1]][2]

features[, `:=`(RESN=paste0(substring(RESN_RESI, 1,1), tolower(substring(RESN_RESI, 2,3))), key_var=paste(uniprot_id, RESI, paste0(substring(RESN_RESI, 1,1), tolower(substring(RESN_RESI, 2,3))), sep=":")  )]

setkey(dbnsfp, key_var)
setkey(features, key_var)

#remove unwanted cols
remove_cols<-colnames(features)[grepl("coord", colnames(features))]
remove_cols<-c(remove_cols, 
               colnames(features)[grepl("RAUTE", colnames(features))])
remove_cols<-c(remove_cols, 
               colnames(features)[grepl("RESN", colnames(features))])
#         -starts_with("coord"),
#         -starts_with("RAUTE"),

merged<-dbnsfp[features, nomatch=0]
merged[,(remove_cols):=NULL]
fwrite(merged, file = args[3])
