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
  mutate(ClinVar2022_pathogenic=(ID %in% clinvar_pathogenic$ID))%>%
  filter(is.na(in_train_ds) | in_train_ds==FALSE)%>%
  filter(!is.na(AlphScore))%>%
  arrange(AlphScore) %>%
  mutate(AlphScore_rank=(1:n())/n())


alphafold_pre_calculated_w_CV2022_train<-alphafold_pre_calculated_w_CV2022 %>%
  filter(in_clinvar_ds==TRUE)%>%
  mutate(outcome=ClinVar2022_pathogenic)

alphafold_pre_calculated_w_CV2022_test<-alphafold_pre_calculated_w_CV2022 %>%
  filter(is.na(in_clinvar_ds) | in_clinvar_ds==FALSE)


set_of_models<-fit_set_of_models(alphafold_pre_calculated_w_CV2022_train)

alphafold_pre_calculated_w_CV2022_test<-predict_set_of_models(set_of_models, alphafold_pre_calculated_w_CV2022_test)


check_count_proteins<-table(alphafold_pre_calculated_w_CV2022_test$Uniprot_acc_split) 
sort(-check_count_proteins)[1:20]

check_count_variants<-table(alphafold_pre_calculated_w_CV2022_test$ID) 
sort(-check_count_variants)[1:20]


roc_rose <- plot(roc(as.integer(alphafold_pre_calculated_w_CV2022_test$ClinVar2022_pathogenic), alphafold_pre_calculated_w_CV2022_test$AlphScore), print.auc = TRUE, 
                 col = "black", print.auc.y = .4)
roc_rose <- plot(roc(alphafold_pre_calculated_w_CV2022_test$ClinVar2022_pathogenic, alphafold_pre_calculated_w_CV2022_test$REVEL_score), print.auc = TRUE, 
                 col = "blue", print.auc.y = .2, add = TRUE)
roc_rose <- plot(roc(alphafold_pre_calculated_w_CV2022_test$ClinVar2022_pathogenic, alphafold_pre_calculated_w_CV2022_test$CADD_raw), print.auc = TRUE, 
                 col = "green", print.auc.y = .6, add = TRUE)
roc_rose <- plot(roc(alphafold_pre_calculated_w_CV2022_test$ClinVar2022_pathogenic, alphafold_pre_calculated_w_CV2022_test$glm_AlphRevel), print.auc = TRUE, 
                 col = "red", print.auc.y = .3, add = TRUE)
roc_rose <- plot(roc(alphafold_pre_calculated_w_CV2022_test$ClinVar2022_pathogenic, alphafold_pre_calculated_w_CV2022_test$glm_AlphCaddDeogen), print.auc = TRUE, 
                 col = "red", print.auc.y = .3, add = TRUE)

ggplot(alphafold_pre_calculated_w_CV2022_test)+
  geom_histogram(aes(x=AlphScore, color=ClinVar2022_pathogenic))


ggplot(alphafold_pre_calculated_w_CV2022_test)+
  geom_point(aes(x=AlphScore, y=CADD_raw, color=ClinVar2022_pathogenic))

