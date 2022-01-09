library(tidyverse)
library(optparse)
library(pROC)

option_list = list(
  make_option(c("-t", "--clinvar_benign"), type="character", default="data/clinvar2022/clinvar_2022_benign.vcf.gz", 
              help="location of ClinVar benign missense variants in vcf like format"),
  make_option(c("-t", "--clinvar_pathogenic"), type="character", default="data/clinvar2022/clinvar_2022_pathogenic.vcf.gz", 
              help="location of ClinVar benign missense variants in vcf like format"),
  make_option(c("-t", "--AlphaFold_scores"), type="character", default="data/merge_all/all_possible_values_concat.csv.gz", 
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
alphafold_pre_calculated<-read_csv(opt$AlphaFold_scores, col_names=) %>%
  mutate(ID=paste(`#chr`, `pos(1-based)`, ref, alt, sep=":"))

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

head(alphafold_pre_calculated)

all_clinvar_ids<-c(clinvar_benign$ID, clinvar_pathogenic$ID)

alphafold_pre_calculated_w_CV2022<-alphafold_pre_calculated %>%
  filter(ID %in% all_clinvar_ids)%>%
  mutate(ClinVar2022_pathogenic=(ID %in% clinvar_pathogenic$ID))%>%
  filter(is.na(in_train_ds) | in_train_ds==FALSE)%>%
  filter(is.na(in_clinvar_ds) | in_clinvar_ds==FALSE)%>%
  filter(!is.na(AlphScore)) %>%
  arrange(AlphScore) %>%
  mutate(AlphScore_rank=(1:n())/n())

alphafold_pre_calculated_w_CV2022

check_count_proteins<-table(alphafold_pre_calculated_w_CV2022$Uniprot_acc_split) 
sort(-check_count_proteins)[1:20]

check_count_variants<-table(alphafold_pre_calculated_w_CV2022$ID) 
sort(-check_count_variants)[1:20]


roc_rose <- plot(roc(as.integer(alphafold_pre_calculated_w_CV2022$ClinVar2022_pathogenic), alphafold_pre_calculated_w_CV2022$AlphScore), print.auc = TRUE, 
                 col = "black", print.auc.y = .4)
roc_rose <- plot(roc(alphafold_pre_calculated_w_CV2022$ClinVar2022_pathogenic, alphafold_pre_calculated_w_CV2022$REVEL_score), print.auc = TRUE, 
                 col = "blue", print.auc.y = .2, add = TRUE)
roc_rose <- plot(roc(alphafold_pre_calculated_w_CV2022$ClinVar2022_pathogenic, alphafold_pre_calculated_w_CV2022$CADD_raw), print.auc = TRUE, 
                 col = "green", print.auc.y = .6, add = TRUE)


ggplot(alphafold_pre_calculated_w_CV2022)+
  geom_histogram(aes(x=AlphScore, color=ClinVar2022_pathogenic))


      