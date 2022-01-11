library(readr)
library(dplyr)
library(optparse)
library(ranger)
source("scripts/existing_scores_glm_functions.R")
set.seed(1)

option_list = list(
  make_option(c("-c", "--csv_location"), type="character", default="data/combine2_protein/P38398_w_AFfeatures.csv.gz", 
              help="csv.gz file"),
  make_option(c("-m", "--model_location"), type="character", default="data/prediction/final_written_full_model.RData", 
              help="Excel file listing columns to use"),
  make_option(c("-o", "--output_file"), type="character", default="data/prediction/predicted.csv.gz", 
              help="Excel file listing columns to use"),
  make_option(c("-u", "--use_cols_file"), type="character", default="data/prediction/final_colnames_to_use.RData", 
              help="Excel file listing columns to use"),
  make_option(c("-t", "--toAS_properties"), type="character", default="data/prediction/final_toAS_properties.RData", 
              help="Excel file listing columns to use"),
  make_option(c("-r", "--reduced"), type="logical", default="FALSE", 
              help="just save essential columns")
)

opt = parse_args(OptionParser(option_list=option_list))

to_AS_table<-read_tsv("resources/to_AS_table.txt")
colnames(to_AS_table) <- paste(colnames(to_AS_table), "toAS", sep = "_")
variants<-read_csv(opt$csv_location, na=c(".","NA", NA))

# if variant file is empty, stop
if (nrow(variants)==0){
  system(paste("touch", opt$output_file))
}else{

model_to_use<-readRDS(opt$model_location)
colnames_to_use<-readRDS(opt$use_cols_file)
toAS_properties<-readRDS(opt$toAS_properties)

variants<-variants %>% 
  left_join(to_AS_table, by=c("to_AS"="to_AS_toAS"))%>%
  mutate(to_AS=toupper(to_AS))%>%
  mutate(from_AS=RESN)%>%
  mutate(var_id_genomic=paste(`#chr`, `pos(1-based)`, sep=":"))%>%
  mutate(var_id_prot=paste(Uniprot_acc_split, RESI, sep=":"))

colnames(variants)<-gsub("++","..",colnames(variants), fixed=TRUE)

variants<-variants%>%
  left_join(toAS_properties, by=c("from_AS"="from_AS_toAS"))

colnames_to_use<-colnames_to_use[colnames_to_use!="outcome"]

sel_vars_to<-colnames(toAS_properties)
sel_vars_to<-sel_vars_to[sel_vars_to!="from_AS_toAS"]
variants[, colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])]<-
  variants[, sel_vars_to] - 
  variants[, colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])]

variants<-variants[complete.cases(variants[,colnames_to_use]),]

variants$AlphScore<-predict(model_to_use, variants)$predictions


variants<-variants %>% 
  mutate(in_train_ds=((is.na(gnomAD_exomes_AC) | (gnomAD_exomes_AC<2 & (gnomAD_exomes_NFE_AC == gnomAD_exomes_AC) ) )  & 
                            gnomAD_genomes_AN>90000 &
                            (!gnomAD_genomes_flag %in% c("lcr","segdup"))&
                            gnomAD_genomes_AC==1 & 
                            gnomAD_genomes_NFE_AC==1) & 
           (is.na(`1000Gp3_AC`) | `1000Gp3_AC`==0)  & 
           (is.na(ESP6500_AA_AC) | ESP6500_AA_AC==0) & 
           (is.na(ESP6500_EA_AC) | ESP6500_AA_AC==0)  |
           gnomAD_genomes_AF> 0.001 & 
           gnomAD_genomes_AN>90000 & 
           (!gnomAD_genomes_flag %in% c("lcr","segdup"))
         )%>%
  mutate(in_clinvar_ds=!is.na(clinvar_clnsig) & clinvar_clnsig %in% c("Likely_pathogenic","Pathogenic","Benign", "Likely_benign") )



if (opt$reduced==TRUE){
 variants<-variants %>% 
  mutate(ID=paste(`#chr`, `pos(1-based)`, ref, alt, sep=":"))%>%
  select(any_of(colnames(variants)[2:10]), 
         ID, genename, Uniprot_acc_split,Uniprot_acc,HGVSp_VEP_split, HGVSp_VEP, CADD_raw, REVEL_score, DEOGEN2_score, 
         b_factor, SOLVENT_ACCESSIBILITY_core, in_train_ds, in_clinvar_ds, AlphScore)
}

write_csv(x=variants,
          file=opt$output_file)
}
