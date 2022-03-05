library(readr)
library(dplyr)
library(optparse)
library(stringr)
source("workflow/scripts/existing_scores_glm_functions.R")
set.seed(1)

option_list = list(
  make_option(c("-t", "--training_dataset"), type="character", default="results/prediction_final/pre_final_model_regular_variants.csv.gz", 
              help="File containing the training dataset"),
  make_option(c("-c", "--combined_model"), type="character", default="results/prediction_final/combined_model.RData", 
              help="File where the final combined logistic regression model should be saved to"),
  make_option(c("-i", "--training_var_ids"), type="character", default="results/prediction_final/training_var_ids.tsv", 
              help="IDs of variants that were in the training / test set")
)

opt = parse_args(OptionParser(option_list=option_list))

variants<-read_csv(opt$training_dataset, na=c(".","NA", NA))
variants$AlphScore<-variants$predicted_Alph

index_col<-get_index_col(variants)

variants$DEOGEN2_score_med<- unlist_score(variants$DEOGEN2_score, index_col)

gnomad_vars<-variants %>% 
  filter(gnomadSet==1)

clinvar_vars_not_gnomad <-variants%>% 
  filter(gnomadSet==0) %>%
  filter(! var_id_genomic %in% gnomad_vars$var_id_genomic )

set_of_models<-fit_set_of_models(clinvar_vars_not_gnomad)

gnomad_set_position<-gnomad_vars$var_id_genomic
clinvar_set_position<-clinvar_vars_not_gnomad$var_id_genomic

training_vars<-rbind(
tibble(var_id_genomic=gnomad_set_position,
       gnomadSet=1),
tibble(var_id_genomic=clinvar_set_position,
       gnomadSet=0)
)
# write to disk
write_tsv(x=training_vars,
          file=opt$training_var_ids)

write_rds(x=set_of_models,
          file=opt$combined_model)

