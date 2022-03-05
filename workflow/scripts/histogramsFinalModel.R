library(tidyverse)
library(optparse)
library(pROC)

source("workflow/scripts/existing_scores_glm_functions.R")
source("workflow/scripts/precision_recall_resource.R")

option_list = list(
  make_option(c("-b", "--clinvar_benign"), type="character", default="results/clinvar2022/clinvar_2022_benign.vcf.gz", 
              help="location of ClinVar benign missense variants in vcf like format"),
  make_option(c("-p", "--clinvar_pathogenic"), type="character", default="results/clinvar2022/clinvar_2022_pathogenic.vcf.gz", 
              help="location of ClinVar benign missense variants in vcf like format"),
  make_option(c("-a", "--AlphaFold_scores"), type="character", default="results/clinvar2022/values_of_clinvar_variants_FINAL.tsv.gz", 
              help="pre calculated Alphafold scores"),
  make_option(c("-v", "--variants"), type="character", default="results/prediction_final/final_regular_variants.csv.gz", 
              help="variants with calculated AlphScore"),
  make_option(c("-o", "--out_folder"), type="character", default="results/final_model_curves/", 
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


all_clinvar_ids<-c(clinvar_benign$ID, clinvar_pathogenic$ID)

alphafold_pre_calculated_w_CV2022<-alphafold_pre_calculated %>%
  filter(ID %in% all_clinvar_ids)%>%
  mutate(outcome=(ID %in% clinvar_pathogenic$ID))%>%
  filter(!is.na(AlphScore))

variants<-variants%>% 
  mutate(ID=paste(`#chr`, `pos(1-based)`,ref, sep=":"))

alphafold_pre_calculated_w_CV2022$pos_in_VEP_and_Uniprot<-get_index_col(alphafold_pre_calculated_w_CV2022)
alphafold_pre_calculated_w_CV2022$DEOGEN2_score_med<-unlist_score(alphafold_pre_calculated_w_CV2022$DEOGEN2_score,
                                                                  alphafold_pre_calculated_w_CV2022$pos_in_VEP_and_Uniprot)

alphafold_pre_calculated_w_CV2022<-alphafold_pre_calculated_w_CV2022 %>%
  mutate(ID=paste(`#chr`, `pos(1-based)`,ref,  sep=":"))%>%
  filter(!ID %in% variants$ID)

# ensure that no variant is multiple times in the data set; if so, take the mean of the predictors, ensure, that the variant is rated consistently benign / pathogenic
testSet<-alphafold_pre_calculated_w_CV2022 %>% 
  group_by(ID) %>%
  summarise(AlphScore=mean(AlphScore), REVEL_score=mean(REVEL_score), glm_AlphRevel=mean(glm_AlphRevel),
            CADD_raw=mean(CADD_raw), glm_AlphCadd=mean(glm_AlphCadd), DEOGEN2_score_med=mean(DEOGEN2_score_med), glm_AlphDeogen=mean(glm_AlphDeogen), outcome=mean(outcome))%>%
  filter(outcome %in% c(0,1))%>%
  mutate(outcome=as.logical(outcome))

nrow(testSet)
sum(complete.cases(testSet))

sum(testSet$ID %in% variants$ID)

testSet<-testSet[complete.cases(testSet),]
table(testSet$outcome)

AlphScorePlot<-ggplot(testSet)+
  geom_density(aes(x=AlphScore, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="AlphScorePlot_FINAL.pdf", 
       plot=AlphScorePlot,
       width=5, height=3)
AlphScore_table<-generate_table("AlphScore")
write_tsv(x=AlphScore_table,
          file="Alphscore_table_FINAL.tsv")

glm_AlphRevel_plot<-ggplot(testSet)+
  geom_density(aes(x=glm_AlphRevel, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="glm_AlphRevel_plot_FINAL.pdf", 
       plot=glm_AlphRevel_plot,
       width=5, height=3)
glm_AlphRevel_table<-generate_table("glm_AlphRevel")
write_tsv(x=glm_AlphRevel_table,
          file="glm_AlphRevel_table_FINAL.tsv")

REVEL_score_plot<-ggplot(testSet)+
  geom_density(aes(x=REVEL_score, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="REVEL_score_plot_FINAL.pdf", 
       plot=REVEL_score_plot,
       width=5, height=3)
REVEL_score_table<-generate_table("REVEL_score")
write_tsv(x=REVEL_score_table,
          file="REVEL_score_table_FINAL.tsv")

CADD_raw_plot<-ggplot(testSet)+
  geom_density(aes(x=CADD_raw, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="CADD_raw_plot_FINAL.pdf", 
       plot=CADD_raw_plot,
       width=6, height=4)
CADD_score_table<-generate_table("CADD_raw")
write_tsv(x=CADD_score_table,
          file="CADD_score_table_FINAL.tsv")

glm_AlphCadd_plot<-ggplot(testSet)+
  geom_density(aes(x=glm_AlphCadd, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="glm_AlphCadd_plot_FINAL.pdf", 
       plot=glm_AlphCadd_plot,
       width=5, height=3)
glm_AlphCadd_table<-generate_table("glm_AlphCadd")
write_tsv(x=glm_AlphCadd_table,
          file="glm_AlphCadd_table_FINAL.tsv")

Deogen_plot<-ggplot(testSet)+
  geom_density(aes(x=DEOGEN2_score_med, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="Deogen_plot_FINAL.pdf", 
       plot=Deogen_plot,
       width=6, height=4)
Deogen_score_table<-generate_table("DEOGEN2_score_med")
write_tsv(x=Deogen_score_table,
          file="Deogen_score_table_FINAL.tsv")

glm_AlphDeogen_plot<-ggplot(testSet)+
  geom_density(aes(x=glm_AlphDeogen, fill=outcome), alpha=0.5)+
  theme_bw()+
  labs(fill = "(likely)\npathogenic:")+
  scale_fill_manual(values=c("darkblue", "red"))
ggsave(filename="glm_AlphDeogen_plot_FINAL.pdf", 
       plot=glm_AlphDeogen_plot,
       width=5, height=3)
glm_AlphDeogen_table<-generate_table("glm_AlphDeogen")
write_tsv(x=glm_AlphDeogen_table,
          file="glm_AlphDeogen_table_FINAL.tsv")



pdf(file="ClinVar_val_REVEL.pdf")
plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
     col = "black", print.auc.y = .35)
plot(roc(testSet$outcome, testSet$REVEL_score), print.auc = TRUE, 
     col = "blue", print.auc.y = .25, add = TRUE)
plot(roc(testSet$outcome, testSet$glm_AlphRevel), print.auc = TRUE, 
     col = "red", print.auc.y = .15, add = TRUE)
legend(0.2, 0.3, legend=c("AlphScore", "REVEL", "AlphScore + REVEL"),
       col=c("black", "blue","red"), lty=1, cex=0.8)
dev.off()

pdf(file="ClinVar_val_CADD.pdf")
plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
     col = "black", print.auc.y = .35)
plot(roc(testSet$outcome, testSet$CADD_raw), print.auc = TRUE, 
     col = "blue", print.auc.y = .25, add = TRUE)
plot(roc(testSet$outcome, testSet$glm_AlphCadd), print.auc = TRUE, 
     col = "red", print.auc.y = .15, add = TRUE)
legend(0.2, 0.3, legend=c("AlphScore", "CADD", "AlphScore + CADD"),
       col=c("black", "blue","red"), lty=1, cex=0.8)
dev.off()

pdf(file="ClinVar_val_DEOGEN2.pdf")
plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
     col = "black", print.auc.y = .35)
plot(roc(testSet$outcome, testSet$DEOGEN2_score_med), print.auc = TRUE, 
     col = "blue", print.auc.y = .25, add = TRUE)
plot(roc(testSet$outcome, testSet$glm_AlphDeogen), print.auc = TRUE, 
     col = "red", print.auc.y = .15, add = TRUE)
legend(0.2, 0.3, legend=c("AlphScore", "DEOGEN2", "AlphScore + DEOGEN2"),
       col=c("black", "blue","red"), lty=1, cex=0.8)
dev.off()

