library(tidyverse)
library(optparse)
library(pROC)
source("scripts/existing_scores_glm_functions.R")

option_list = list(
  make_option(c("-t", "--clinvar_benign"), type="character", default="data/clinvar2022/clinvar_2022_benign.vcf.gz", 
              help="location of ClinVar benign missense variants in vcf like format"),
  make_option(c("-t", "--clinvar_pathogenic"), type="character", default="data/clinvar2022/clinvar_2022_pathogenic.vcf.gz", 
              help="location of ClinVar benign missense variants in vcf like format"),
  make_option(c("-t", "--AlphaFold_scores"), type="character", default="data/clinvar2022/values_of_clinvar_variants.tsv.gz", 
              help="pre calculated Alphafold scores"),
  make_option(c("-t", "--header_alph_scores"), type="character", default="data/merge_all/header.csv", 
              help="pre calculated Alphafold scores"),
  make_option(c("-o", "--out_folder"), type="character", default="data/clinvar2022_alphafold/", 
              help="name of folder to store output")
              )
opt = parse_args(OptionParser(option_list=option_list))

# load data, add variant ID (chrom:pos:ref:alt)
clinvar_benign<-read_tsv(opt$clinvar_benign, col_names=TRUE) %>%
  mutate(ID=paste(CHROM,POS,REF,ALT, sep=":"))
clinvar_pathogenic<-read_tsv(opt$clinvar_pathogenic, col_names=TRUE) %>%
  mutate(ID=paste(CHROM,POS,REF,ALT, sep=":"))
header_alph<-read_csv(opt$header_alph_scores)
alphafold_pre_calculated<-read_tsv(opt$AlphaFold_scores, col_names=colnames(header_alph))%>%
  mutate(ID=paste(`#chr`, `pos(1-based)`, ref, alt, sep=":"))

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

head(alphafold_pre_calculated)
alphafold_pre_calculated

alphafold_pre_calculated$pos_in_VEP_and_Uniprot<-get_index_col(alphafold_pre_calculated)
alphafold_pre_calculated$DEOGEN2_score_med<-unlist_score(alphafold_pre_calculated$DEOGEN2_score, alphafold_pre_calculated$pos_in_VEP_and_Uniprot)


all_clinvar_ids<-c(clinvar_benign$ID, clinvar_pathogenic$ID)

alphafold_pre_calculated_w_CV2022<-alphafold_pre_calculated %>%
  filter(ID %in% all_clinvar_ids)%>%
  mutate(outcome=(ID %in% clinvar_pathogenic$ID))%>%
  filter(!is.na(AlphScore))

num_vars<-nrow(alphafold_pre_calculated_w_CV2022)
alphafold_pre_calculated_w_CV2022$CrossValGroup<-as.integer(runif(num_vars,1,6))


score_performance_tbl<-tibble()

for (i in unique(alphafold_pre_calculated_w_CV2022$CrossValGroup)){
  
  trainSet<-alphafold_pre_calculated_w_CV2022 %>%
    filter(CrossValGroup!=i)
  
  testSet<-alphafold_pre_calculated_w_CV2022 %>%
    filter(CrossValGroup==i)%>%
    filter(!ID %in% trainSet$ID)
  
  set_of_models<-fit_set_of_models(trainSet)
  testSet<-predict_set_of_models(set_of_models, testSet)
  
  roc_rose <- plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
                   col = "black", print.auc.y = .4)
  roc_rose <- plot(roc(testSet$outcome, testSet$REVEL_score), print.auc = TRUE, 
                   col = "blue", print.auc.y = .2, add = TRUE)
  roc_rose <- plot(roc(testSet$outcome, testSet$glm_AlphRevel), print.auc = TRUE, 
                   col = "red", print.auc.y = .3, add = TRUE)
  roc_rose <- plot(roc(testSet$outcome, testSet$Alph_null), print.auc = TRUE, 
                   col = "grey", print.auc.y = .1, add = TRUE)
  
  roc_rose <- plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
                   col = "black", print.auc.y = .4)
  roc_rose <- plot(roc(testSet$outcome, testSet$CADD_raw), print.auc = TRUE, 
                   col = "blue", print.auc.y = .2, add = TRUE)
  roc_rose <- plot(roc(testSet$outcome, testSet$glm_AlphCadd), print.auc = TRUE, 
                   col = "green", print.auc.y = .6, add = TRUE)
  roc_rose <- plot(roc(testSet$outcome, testSet$Alph_null), print.auc = TRUE, 
                   col = "grey", print.auc.y = .1, add = TRUE)
  
  roc_rose <- plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
                   col = "black", print.auc.y = .4)
  roc_rose <- plot(roc(testSet$outcome, testSet$DEOGEN2_score_med), print.auc = TRUE, 
                   col = "blue", print.auc.y = .2, add = TRUE)
  roc_rose <- plot(roc(testSet$outcome, testSet$glm_AlphDeogen), print.auc = TRUE, 
                   col = "green", print.auc.y = .6, add = TRUE)
  roc_rose <- plot(roc(testSet$outcome, testSet$Alph_null), print.auc = TRUE, 
                   col = "grey", print.auc.y = .1, add = TRUE)
  
  score_performance_tbl <- rbind(score_performance_tbl, 
                       tibble(auc_Alph=roc(testSet$outcome, testSet$AlphScore)$auc, 
                              auc_Alph_null=roc(testSet$outcome, testSet$Alph_null)$auc, 
                              auc_CADD=roc(testSet$outcome, testSet$CADD_raw)$auc, 
                              auc_REVEL=roc(testSet$outcome, testSet$REVEL_score)$auc,
                              auc_DEOGEN2=roc(testSet$outcome, testSet$DEOGEN2_score_med)$auc, 
                              auc_AlphCadd=roc(testSet$outcome, testSet$glm_AlphCadd)$auc, 
                              auc_AlphDeogen=roc(testSet$outcome, testSet$glm_AlphDeogen)$auc, 
                              auc_AlphRevel=roc(testSet$outcome, testSet$glm_AlphRevel)$auc, 
                              num_train=nrow(trainSet),
                              num_test=nrow(testSet) ))
}

score_mean<-score_performance_tbl %>%
  summarise_all(mean)

score_sd<-score_performance_tbl %>%
  summarise_all(sd)

return_sem<-function(x) { return(sd(x)/sqrt(length(x)) )}
score_se<-score_performance_tbl %>%
  summarise_all(return_sem)

