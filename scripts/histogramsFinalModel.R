library(tidyverse)
library(optparse)

source("scripts/existing_scores_glm_functions.R")
source("scripts/precision_recall_resource.R")

option_list = list(
  make_option(c("-b", "--clinvar_benign"), type="character", default="data/clinvar2022/clinvar_2022_benign.vcf.gz", 
              help="location of ClinVar benign missense variants in vcf like format"),
  make_option(c("-p", "--clinvar_pathogenic"), type="character", default="data/clinvar2022/clinvar_2022_pathogenic.vcf.gz", 
              help="location of ClinVar benign missense variants in vcf like format"),
  make_option(c("-a", "--AlphaFold_scores"), type="character", default="data/clinvar2022/values_of_clinvar_variants.tsv.gz", 
              help="pre calculated Alphafold scores"),
  make_option(c("-v", "--variants"), type="character", default="data/prediction_final/pre_final_model_regular_variants.csv.gz", 
              help="variants with calculated AlphScore"),
  make_option(c("-o", "--out_folder"), type="character", default="data/clinvar2022_alphafold/", 
              help="name of folder to store output")
              )
opt = parse_args(OptionParser(option_list=option_list))

# load data, add variant ID (chrom:pos:ref:alt)
clinvar_benign<-read_tsv(opt$clinvar_benign, col_names=TRUE) %>%
  mutate(ID=paste(CHROM,POS,REF,ALT, sep=":"))
clinvar_pathogenic<-read_tsv(opt$clinvar_pathogenic, col_names=TRUE) %>%
  mutate(ID=paste(CHROM,POS,REF,ALT, sep=":"))
alphafold_pre_calculated<-read_tsv(opt$AlphaFold_scores, col_names=TRUE, 
                                   col_types = cols(.default = "?", REVEL_score = "d"))%>%
  mutate(ID=paste(`#chr`, `pos(1-based)`, ref, alt, sep=":"))
variants<-read_csv(opt$variants, 
                   col_types = cols(.default = "?", REVEL_score = "d"))


dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

head(alphafold_pre_calculated)

alphafold_pre_calculated$pos_in_VEP_and_Uniprot<-get_index_col(alphafold_pre_calculated)
alphafold_pre_calculated$DEOGEN2_score_med<-unlist_score(alphafold_pre_calculated$DEOGEN2_score, alphafold_pre_calculated$pos_in_VEP_and_Uniprot)


all_clinvar_ids<-c(clinvar_benign$ID, clinvar_pathogenic$ID)

alphafold_pre_calculated_w_CV2022<-alphafold_pre_calculated %>%
  filter(ID %in% all_clinvar_ids)%>%
  mutate(outcome=(ID %in% clinvar_pathogenic$ID))%>%
  filter(!is.na(AlphScore))


variants$AlphScore<-variants$predicted_Alph


variants<-variants%>% 
  mutate(ID=paste(`#chr`, `pos(1-based)`,ref, sep=":"))

### combine scores on interim dataset:
interim_dataset<-variants %>% 
  filter(clinvar_interim_test==TRUE)
gnomad_dataset<-variants %>% 
  filter(gnomad_train==TRUE)
non_test_dataset<-variants %>% 
  filter(clinvar_holdout_test==FALSE)
test_dataset<-variants %>% 
  filter(clinvar_holdout_test==TRUE)

alphafold_pre_calculated_w_CV2022<-alphafold_pre_calculated_w_CV2022 %>%
  mutate(ID=paste(`#chr`, `pos(1-based)`,ref,  sep=":"))%>%
  filter(!ID %in% interim_dataset$ID)%>%
  filter(!ID %in% non_test_dataset$ID)

set_of_models<-fit_set_of_models(interim_dataset)
testSet<-predict_set_of_models(set_of_models, alphafold_pre_calculated_w_CV2022)

# ensure that no variant is multiple times in the data set; if so, take the mean of the predictors, ensure, that the variant is rated consistently benign / pathogenic
testSet<-testSet %>% 
  group_by(ID) %>%
  summarise(AlphScore=mean(AlphScore), REVEL_score=mean(REVEL_score), glm_AlphRevel=mean(glm_AlphRevel), Alph_null=mean(Alph_null),
            CADD_raw=mean(CADD_raw), glm_AlphCadd=mean(glm_AlphCadd), DEOGEN2_score_med=mean(DEOGEN2_score_med), glm_AlphDeogen=mean(glm_AlphDeogen), outcome=mean(outcome))%>%
  filter(outcome %in% c(0,1))%>%
  mutate(outcome=as.logical(outcome))

nrow(testSet)
sum(complete.cases(testSet))

sum(testSet$ID %in% test_dataset$ID)
sum(testSet$ID %in% variants$ID)

testSet<-testSet[complete.cases(testSet),]
table(testSet$outcome)

AlphNullPlot<-ggplot(testSet)+
  geom_density(aes(x=Alph_null, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="AlphNullPlot.pdf", 
       plot=AlphNullPlot,
       width=5, height=3)
alph_null_table<-generate_table("Alph_null")
write_tsv(x=alph_null_table,
          file="alph_null_table.tsv")

AlphScorePlot<-ggplot(testSet)+
  geom_density(aes(x=AlphScore, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="AlphScorePlot.pdf", 
       plot=AlphScorePlot,
       width=5, height=3)
AlphScore_table<-generate_table("AlphScore")
write_tsv(x=alph_null_table,
          file="alph_null_table.tsv")

glm_AlphRevel_plot<-ggplot(testSet)+
  geom_density(aes(x=glm_AlphRevel, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="glm_AlphRevel_plot.pdf", 
       plot=glm_AlphRevel_plot,
       width=5, height=3)
glm_AlphRevel_table<-generate_table("glm_AlphRevel")
write_tsv(x=glm_AlphRevel_table,
          file="glm_AlphRevel_table.tsv")

REVEL_score_plot<-ggplot(testSet)+
  geom_density(aes(x=REVEL_score, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="REVEL_score_plot.pdf", 
       plot=REVEL_score_plot,
       width=5, height=3)
REVEL_score_table<-generate_table("REVEL_score")
write_tsv(x=REVEL_score_table,
          file="REVEL_score_table.tsv")

CADD_raw_plot<-ggplot(testSet)+
  geom_density(aes(x=CADD_raw, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="CADD_raw_plot.pdf", 
       plot=CADD_raw_plot,
       width=6, height=4)
CADD_score_table<-generate_table("CADD_raw")
write_tsv(x=CADD_score_table,
          file="CADD_score_table.tsv")

glm_AlphCadd_plot<-ggplot(testSet)+
  geom_density(aes(x=glm_AlphCadd, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="glm_AlphCadd_plot.pdf", 
       plot=glm_AlphCadd_plot,
       width=5, height=3)
glm_AlphCadd_table<-generate_table("glm_AlphCadd")
write_tsv(x=glm_AlphCadd_table,
          file="glm_AlphCadd_table.tsv")
