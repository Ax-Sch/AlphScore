library(tidyverse)
library(optparse)
library(pROC)
library(PRROC)
library(VennDiagram)
library(boot)

source("scripts/existing_scores_glm_functions.R")
source("scripts/precision_recall_resource.R")

BOOT_REPETITIONS=1000
set.seed(1)

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

variants$pos_in_VEP_and_Uniprot<-get_index_col(variants)
variants$DEOGEN2_score_med<-unlist_score(variants$DEOGEN2_score, variants$pos_in_VEP_and_Uniprot)

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



score_performance_tbl<-tibble()

pdf(file="ClinVar_val_REVEL.pdf")
plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
       col = "black", print.auc.y = .35)
  plot(roc(testSet$outcome, testSet$REVEL_score), print.auc = TRUE, 
       col = "blue", print.auc.y = .25, add = TRUE)
  plot(roc(testSet$outcome, testSet$glm_AlphRevel), print.auc = TRUE, 
       col = "red", print.auc.y = .15, add = TRUE)
  plot(roc(testSet$outcome, testSet$Alph_null), print.auc = TRUE, 
       col = "grey", print.auc.y = .05, add = TRUE)
  legend(0.2, 0.3, legend=c("AlphScore", "REVEL", "AlphScore + REVEL", "NullModel"),
         col=c("black", "blue","red", "grey"), lty=1, cex=0.8)
dev.off()

pdf(file="ClinVar_val_CADD.pdf")
  plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
       col = "black", print.auc.y = .35)
  plot(roc(testSet$outcome, testSet$CADD_raw), print.auc = TRUE, 
       col = "blue", print.auc.y = .25, add = TRUE)
  plot(roc(testSet$outcome, testSet$glm_AlphCadd), print.auc = TRUE, 
       col = "red", print.auc.y = .15, add = TRUE)
  plot(roc(testSet$outcome, testSet$Alph_null), print.auc = TRUE, 
       col = "grey", print.auc.y = .05, add = TRUE)
  legend(0.2, 0.3, legend=c("AlphScore", "CADD", "AlphScore + CADD", "NullModel"),
         col=c("black", "blue","red", "grey"), lty=1, cex=0.8)
dev.off()
  
pdf(file="ClinVar_val_DEOGEN2.pdf")
  plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
       col = "black", print.auc.y = .35)
  plot(roc(testSet$outcome, testSet$DEOGEN2_score_med), print.auc = TRUE, 
       col = "blue", print.auc.y = .25, add = TRUE)
  plot(roc(testSet$outcome, testSet$glm_AlphDeogen), print.auc = TRUE, 
       col = "red", print.auc.y = .15, add = TRUE)
  plot(roc(testSet$outcome, testSet$Alph_null), print.auc = TRUE, 
       col = "grey", print.auc.y = .05, add = TRUE)
  legend(0.2, 0.3, legend=c("AlphScore", "DEOGEN2", "AlphScore + DEOGEN2", "NullModel"),
         col=c("black", "blue","red", "grey"), lty=1, cex=0.8)
dev.off()
  




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


pr_cadd<-pr.curve(scores.class0=testSet$CADD_raw, weights.class0=testSet$outcome, curve=T)
pr_alphCadd<-pr.curve(scores.class0=testSet$glm_AlphCadd, weights.class0=testSet$outcome, curve=T)
pr_alph<-pr.curve(scores.class0=testSet$AlphScore, weights.class0=testSet$outcome, curve=T)
pr_alph_null<-pr.curve(scores.class0=testSet$Alph_null, weights.class0=testSet$outcome, curve=T)

pdf(file="proc_auc_CADD.pdf")
  plot(pr_alphCadd, color="red")
  plot(pr_cadd, add=TRUE, color="blue")
  plot(pr_alph, add=TRUE, color="black")
  plot(pr_alph_null, add=TRUE, color="grey")
  legend(0.2, 0.3, legend=c("AlphScore", "CADD", "AlphScore + CADD", "NullModel"),
         col=c("black", "blue","red", "grey"), lty=1, cex=0.8)
dev.off()


pr_revel<-pr.curve(scores.class0=testSet$REVEL_score, weights.class0=testSet$outcome, curve=T)
pr_alphRevel<-pr.curve(scores.class0=testSet$glm_AlphRevel, weights.class0=testSet$outcome, curve=T)

pdf(file="proc_auc_REVEL.pdf")
  plot(pr_alphRevel, color="red")
  plot(pr_revel, add=TRUE, color="blue")
  plot(pr_alph, add=TRUE, color="black")
  plot(pr_alph_null, add=TRUE, color="grey")
  legend(0.2, 0.3, legend=c("AlphScore", "REVEL", "AlphScore + REVEL", "NullModel"),
         col=c("black", "blue","red", "grey"), lty=1, cex=0.8)
dev.off()

pr_deogen<-pr.curve(scores.class0=testSet$DEOGEN2_score_med, weights.class0=testSet$outcome, curve=T)
pr_alphDeogen<-pr.curve(scores.class0=testSet$glm_AlphDeogen, weights.class0=testSet$outcome, curve=T)

pdf(file="proc_auc_DEOGEN2.pdf")
  plot(pr_alphDeogen, color="red")
  plot(pr_deogen, add=TRUE, color="blue")
  plot(pr_alph, add=TRUE, color="black")
  plot(pr_alph_null, add=TRUE, color="grey")
  legend(0.2, 0.3, legend=c("AlphScore", "DEOGEN2", "AlphScore + DEOGEN2", "NullModel"),
         col=c("black", "blue","red", "grey"), lty=1, cex=0.8)
dev.off()



# get aucs and do bootstrapping
get_score_performance_table<-function(testSet_private, i){
  testSet_private<-testSet_private[i,]
  score_performance_tbl_private <-
    tibble(
      Alph_ROC=roc(testSet_private$outcome, testSet_private$AlphScore)$auc, 
      Alph_null_ROC=roc(testSet_private$outcome, testSet_private$Alph_null)$auc, 
      CADD_ROC=roc(testSet_private$outcome, testSet_private$CADD_raw)$auc, 
      REVEL_ROC=roc(testSet_private$outcome, testSet_private$REVEL_score)$auc,
      DEOGEN2_ROC=roc(testSet_private$outcome, testSet_private$DEOGEN2_score_med)$auc, 
      AlphCadd_ROC=roc(testSet_private$outcome, testSet_private$glm_AlphCadd)$auc, 
      AlphDeogen_ROC=roc(testSet_private$outcome, testSet_private$glm_AlphDeogen)$auc, 
      AlphRevel_ROC=roc(testSet_private$outcome, testSet_private$glm_AlphRevel)$auc, 
      
      Alph_PROC=pr.curve(weights.class0=testSet_private$outcome, scores.class0=testSet_private$AlphScore)$auc.integral, 
      Alph_null_PROC=pr.curve(weights.class0=testSet_private$outcome, scores.class0=testSet_private$Alph_null)$auc.integral, 
      CADD_PROC=pr.curve(weights.class0=testSet_private$outcome, scores.class0=testSet_private$CADD_raw)$auc.integral, 
      REVEL_PROC=pr.curve(weights.class0=testSet_private$outcome, scores.class0=testSet_private$REVEL_score)$auc.integral,
      DEOGEN2_PROC=pr.curve(weights.class0=testSet_private$outcome, scores.class0=testSet_private$DEOGEN2_score_med)$auc.integral, 
      AlphCadd_PROC=pr.curve(weights.class0=testSet_private$outcome, scores.class0=testSet_private$glm_AlphCadd)$auc.integral, 
      AlphDeogen_PROC=pr.curve(weights.class0=testSet_private$outcome, scores.class0=testSet_private$glm_AlphDeogen)$auc.integral, 
      AlphRevel_PROC=pr.curve(weights.class0=testSet_private$outcome, scores.class0=testSet_private$glm_AlphRevel)$auc.integral,             
      
      num_test=nrow(testSet_private) )
  
  return(score_performance_tbl_private)
}


get_score_performance_table_boot<-function(testSet_private,i){
  return(unlist(get_score_performance_table(testSet_private,i)))
}

score_performance_tbl<-get_score_performance_table(testSet, TRUE)

write_tsv(x=score_performance_tbl, 
          file="score_performance_tbl.tsv")


compare_scores<-function(booted_values_private){
  booted_values_private_summarised<-booted_values_private %>%
    rowwise%>%
    mutate(rowmax=max(across()))%>%
    mutate(AlphCadd_biggerThanCadd_ROC=AlphCadd_ROC>CADD_ROC)%>%
    mutate(AlphDeogen_biggerThanDeogen_ROC=AlphDeogen_ROC>DEOGEN2_ROC)%>%
    mutate(AlphRevel_biggerThanRevel_ROC=AlphRevel_ROC>REVEL_ROC)%>%
    mutate(AlphCadd_biggerThanCadd_PROC=AlphCadd_PROC>CADD_PROC)%>%
    mutate(AlphDeogen_biggerThanDeogen_PROC=AlphDeogen_PROC>DEOGEN2_PROC)%>%
    mutate(AlphRevel_biggerThanRevel_PROC=AlphRevel_PROC>REVEL_PROC)
  return(booted_values_private_summarised)
}

get_p_value_table<-function(booted_values_summarised_private){
  relevant_scores<-c("AlphDeogen_biggerThanDeogen_ROC","AlphRevel_biggerThanRevel_ROC","AlphCadd_biggerThanCadd_ROC",
                     "AlphDeogen_biggerThanDeogen_PROC","AlphRevel_biggerThanRevel_PROC","AlphCadd_biggerThanCadd_PROC")
  p_val_tabl<-tibble()
  for (columname in relevant_scores){
    p_val_tabl<- rbind(p_val_tabl, 
                       tibble(name=columname, num_true=sum(unlist(booted_values_summarised_private[,columname])), total_num=nrow(booted_values_summarised_private)))
  }
  return(p_val_tabl)
}

boot_aucs<-boot(data=testSet, statistic=get_score_performance_table_boot, R=BOOT_REPETITIONS, parallel="multicore", ncpus=3)
booted_aucs_table<-as_tibble(boot_aucs$t)
colnames(booted_aucs_table)<-colnames(score_performance_tbl)
booted_aucs_table<-compare_scores(booted_aucs_table)
pValTable_AbsCor<-get_p_value_table(booted_aucs_table)
pValTable_AbsCor<-pValTable_AbsCor%>% 
  mutate(pval=1-num_true/total_num)
pValTable_AbsCor

write_tsv(x = pValTable_AbsCor,
          file = "pValTable_AbsCor.tsv")



# Diagram with ROCs

score_performance_tbl_spread<-score_performance_tbl%>%
  select(-num_test)%>%
  gather(key="method")%>%
  arrange(value)

score_performance_tbl_spread<-score_performance_tbl_spread%>%
  mutate(method=factor(method, levels=c("Alph_null_ROC","Alph_ROC","CADD_ROC","AlphCadd_ROC","DEOGEN2_ROC","AlphDeogen_ROC","REVEL_ROC","AlphRevel_ROC")))

plot_aucs_ClinVar<-ggplot(score_performance_tbl_spread, aes(x=method, y=value))+
  stat_summary(fun.y = mean, geom = "bar") + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=10),
        axis.text.y = element_text(color="black", size=10))+
  coord_cartesian(ylim=c(0.5,1)) + 
  labs(x = "")+
  labs(y = "AUC (ClinVar test set)", size=12)

plot_aucs_ClinVar
ggsave(filename= "plot_aucs_ClinVar.pdf", plot=plot_aucs_ClinVar, height=5, width=4)







datasets_venn<-venn.diagram(x=list(gnomad_dataset$ID, interim_dataset$ID, testSet$ID), 
             category.names = c("gnomAD_train" , "ClinVar_validation", "ClinVar_test/hold-out"), 
             filename = NULL)

ggsave(datasets_venn, file="datasets_venn.pdf", device = "pdf", width=6, height=6)

