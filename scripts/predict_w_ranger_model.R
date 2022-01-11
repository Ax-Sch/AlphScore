library(readr)
library(dplyr)
library(stringr)
library(optparse)
library(ranger)
source("scripts/existing_scores_glm_functions.R")
set.seed(1)

option_list = list(
  make_option(c("-c", "--csv_location"), type="character", default="data/combine2_protein/P38398_w_AFfeatures.csv.gz", 
              help="csv.gz file"),
  make_option(c("-m", "--model_location"), type="character", default="data/prediction/final_written_full_model.RData", 
              help="File containing the full model for prediction"),
  make_option(c("-u", "--use_cols_file"), type="character", default="data/prediction/final_colnames_to_use.RData", 
              help="File containing the colnames to use in the final model"),
  make_option(c("-t", "--toAS_properties"), type="character", default="data/prediction/final_toAS_properties.RData", 
              help="File containing the properties of the alternative amino acids"),
  make_option(c("-n", "--model_location_null"), type="character", default="data/prediction/final_written_full_model.RData", 
              help="File containing the full model for prediction, null model"),
  make_option(c("-x", "--use_cols_file_null"), type="character", default="data/prediction/final_colnames_to_use.RData", 
              help="File containing the colnames to use in the final model, null model"),
  make_option(c("-v", "--toAS_properties_null"), type="character", default="data/prediction/final_toAS_properties.RData", 
              help="File containing the properties of the alternative amino acids, null model"),
  make_option(c("-r", "--reduced"), type="logical", default="FALSE", 
              help="just save essential columns"),
  make_option(c("-o", "--output_file"), type="character", default="data/prediction/predicted.csv.gz", 
              help="Excel file listing columns to use")
)

opt = parse_args(OptionParser(option_list=option_list))

to_AS_table<-read_tsv("resources/to_AS_table.txt")
variants<-read_csv(opt$csv_location, na=c(".","NA", NA))

# if variant file is empty, stop
if (nrow(variants)==0){
  system(paste("touch", opt$output_file))
}else{

model_to_use<-readRDS(opt$model_location)
colnames_to_use<-readRDS(opt$use_cols_file)
toAS_properties<-readRDS(opt$toAS_properties)

model_to_use_null<-readRDS(opt$model_location_null)
colnames_to_use_null<-readRDS(opt$use_cols_file_null)
toAS_properties_null<-readRDS(opt$toAS_properties_null)


# function to add to AS properties
addToAS<-function(variants_par, toAsProp, colToUse){
variants_mod<-variants
variants_mod<-variants_mod%>%
    left_join(toAsProp, by=c("from_AS"="from_AS_toAS"))
  
  sel_vars_to<-colnames(toAsProp)
  sel_vars_to<-sel_vars_to[sel_vars_to!="from_AS_toAS"]
  sel_vars_from<-str_replace(sel_vars_to, fixed("_toAS"),"")
  
  variants_mod[, sel_vars_from]<-variants_mod[, sel_vars_to] - variants_mod[, sel_vars_from]
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

# predict null model
variants_alph_null<-addToAS(variants, toAS_properties_null, colnames_to_use_null)
variants$Alph_null<-predict(model_to_use_null, variants_alph_null)$predictions

# flag variants that are in the training data set or that are in dbNSFP ClinVar
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


# if wanted, reduce the output to have smaller files
if (opt$reduced==TRUE){
 variants<-variants %>% 
  mutate(ID=paste(`#chr`, `pos(1-based)`, ref, alt, sep=":"))%>%
  select(any_of(colnames(variants)[2:10]), 
         ID, genename, Uniprot_acc_split,Uniprot_acc,HGVSp_VEP_split, HGVSp_VEP, CADD_raw, REVEL_score, DEOGEN2_score, 
         b_factor, SOLVENT_ACCESSIBILITY_core, in_train_ds, in_clinvar_ds, Alph_null, AlphScore)
}

# write to disk
write_csv(x=variants,
          file=opt$output_file)
}

