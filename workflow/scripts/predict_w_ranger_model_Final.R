library(readr)
library(dplyr)
library(stringr)
library(optparse)
library(ranger)
source("workflow/scripts/existing_scores_glm_functions.R")
set.seed(1)

option_list = list(
  make_option(c("-c", "--csv_location"), type="character", default="results/combine2_protein/A0A1B0GUS4_w_AFfeatures.csv.gz", 
              help="csv.gz file"),
  make_option(c("-m", "--model_location"), type="character", default="results/prediction/final_written_full_model.RData", 
              help="File containing the full model for prediction"),
  make_option(c("-u", "--use_cols_file"), type="character", default="results/prediction/final_colnames_to_use.RData", 
              help="File containing the colnames to use in the final model"),
  make_option(c("-t", "--toAS_properties"), type="character", default="results/prediction/final_toAS_properties.RData", 
              help="File containing the properties of the alternative amino acids"),
  make_option(c("-b", "--combined_model"), type="character", default="results/prediction_final/combined_model.RData", 
              help="File where the final combined logistic regression model should be saved to"),
  make_option(c("-v", "--training_var_ids"), type="character", default="results/prediction_final/training_var_ids.tsv", 
              help="IDs of variants that were in the training / test set"),
  make_option(c("-o", "--output_file"), type="character", default="results/prediction/predicted.csv.gz", 
              help="Excel file listing columns to use"),
  make_option(c("-r", "--reduced"), type="logical", default="FALSE", 
              help="just save essential columns")
)

opt = parse_args(OptionParser(option_list=option_list))

to_AS_table<-read_tsv("resources/to_AS_table.txt")
variants<-read_csv(opt$csv_location, na=c(".","NA", NA), col_types = cols(.default = "?", REVEL_score = "d"))

# if variant file is empty, stop
if (nrow(variants)==0){
  system(paste("touch", opt$output_file))
}else{

model_to_use<-readRDS(opt$model_location)
colnames_to_use<-readRDS(opt$use_cols_file)
toAS_properties<-readRDS(opt$toAS_properties)

set_of_models<-readRDS(opt$combined_model)

training_var_ids<-read_tsv(opt$training_var_ids)
gnomad_var_ids<-(training_var_ids %>% filter(gnomadSet==1))$var_id_genomic
ClinVar_var_ids<-(training_var_ids %>% filter(gnomadSet==0))$var_id_genomic

# function to add to AS properties
addToAS<-function(variants_par, toAsProp, colToUse){
  variants_mod<-variants
  variants_mod<-variants_mod%>%
    left_join(toAsProp, by=c("from_AS"="from_AS_toAS"))
  
  colnames_toAS<-colnames(toAsProp)
  colnames_toAS<-colnames_toAS[colnames_toAS!="from_AS_toAS"]
  
  sel_vars_to<-str_replace(colnames_toAS, fixed("_toAS"),"")
  
  variants_mod[, colnames_toAS]<-variants_mod[, sel_vars_to] - variants_mod[, colnames_toAS]
  colToUse_mod<-colToUse[colToUse!="outcome"]
  
  return(variants_mod[, colToUse_mod])
}


# external function that is also used by preprocess
variants<-prepareVariantsForPrediction(variants, to_AS_table)

#add to AS properties
variants_alph<-addToAS(variants, toAS_properties, colnames_to_use)

# just keep complete cases
variants<-variants[complete.cases(variants_alph),]
variants_alph<-variants_alph[complete.cases(variants_alph),]

# predict AlphScore
variants$AlphScore<-predict(model_to_use, variants_alph)$predictions
rm(variants_alph)

index_col<-get_index_col(variants)

variants$DEOGEN2_score_med<- unlist_score(variants$DEOGEN2_score, index_col)

variants<-predict_set_of_models(set_of_models = set_of_models, variants_to_predict = variants)
colNamesCombinedModels<-unlist(set_of_models[2])

# flag variants that are in the training data set or that are in dbNSFP ClinVar
variants<-variants %>% 
  mutate(in_gnomad_train=var_id_genomic %in% gnomad_var_ids)%>%
  mutate(in_clinvar_ds=var_id_genomic %in% ClinVar_var_ids)


# if wanted, reduce the output to have smaller files
if (opt$reduced==TRUE){
 variants<-variants %>% 
  mutate(ID=paste(`#chr`, `pos(1-based)`, ref, alt, sep=":"))%>%
  select(any_of(colnames(variants)[2:10]), 
         ID, genename, Uniprot_acc_split,Uniprot_acc,HGVSp_VEP_split, HGVSp_VEP, CADD_raw, REVEL_score, DEOGEN2_score, 
         b_factor, SOLVENT_ACCESSIBILITY_core, in_gnomad_train, in_clinvar_ds, AlphScore, any_of(colNamesCombinedModels))
}

# write to disk
write_csv(x=variants,
          file=opt$output_file)
}

