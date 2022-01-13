library(tidyverse)
library(optparse)
library(pROC)
source("scripts/existing_scores_glm_functions.R")

option_list = list(
  make_option(c("-b", "--clinvar_benign"), type="character", default="data/clinvar2022/clinvar_2022_benign.vcf.gz", 
              help="location of ClinVar benign missense variants in vcf like format"),
  make_option(c("-p", "--clinvar_pathogenic"), type="character", default="data/clinvar2022/clinvar_2022_pathogenic.vcf.gz", 
              help="location of ClinVar benign missense variants in vcf like format"),
  make_option(c("-a", "--AlphaFold_scores"), type="character", default="data/clinvar2022/values_of_clinvar_variants.tsv.gz", 
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
alphafold_pre_calculated<-read_tsv(opt$AlphaFold_scores, col_names=TRUE)%>%
  mutate(ID=paste(`#chr`, `pos(1-based)`, ref, alt, sep=":"))

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

head(alphafold_pre_calculated)

alphafold_pre_calculated$pos_in_VEP_and_Uniprot<-get_index_col(alphafold_pre_calculated)
alphafold_pre_calculated$DEOGEN2_score_med<-unlist_score(alphafold_pre_calculated$DEOGEN2_score, alphafold_pre_calculated$pos_in_VEP_and_Uniprot)

all_clinvar_ids<-c(clinvar_benign$ID, clinvar_pathogenic$ID)

alphafold_pre_calculated_w_CV2022<-alphafold_pre_calculated %>%
  filter(!in_gnomad_train | is.na(in_gnomad_train))%>%
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
  
  plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
       col = "black", print.auc.y = .35)
  plot(roc(testSet$outcome, testSet$REVEL_score), print.auc = TRUE, 
       col = "blue", print.auc.y = .25, add = TRUE)
  plot(roc(testSet$outcome, testSet$glm_AlphRevel), print.auc = TRUE, 
       col = "red", print.auc.y = .15, add = TRUE)
  plot(roc(testSet$outcome, testSet$Alph_null), print.auc = TRUE, 
       col = "grey", print.auc.y = .05, add = TRUE)
  legend(0.2, 0.3, legend=c("AlphScore", "REVEL", "AlphRevel", "Alph_null"),
         col=c("black", "blue","red", "grey"), lty=1, cex=0.8)
  title(main = paste("Hold-out-set number", i))
  
  plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
       col = "black", print.auc.y = .35)
  plot(roc(testSet$outcome, testSet$CADD_raw), print.auc = TRUE, 
       col = "blue", print.auc.y = .25, add = TRUE)
  plot(roc(testSet$outcome, testSet$glm_AlphCadd), print.auc = TRUE, 
       col = "red", print.auc.y = .15, add = TRUE)
  plot(roc(testSet$outcome, testSet$Alph_null), print.auc = TRUE, 
       col = "grey", print.auc.y = .05, add = TRUE)
  legend(0.2, 0.3, legend=c("AlphScore", "CADD", "AlphCadd", "Alph_null"),
         col=c("black", "blue","red", "grey"), lty=1, cex=0.8)
  title(main = paste("Hold-out-set number", i))
  
  plot(roc(as.integer(testSet$outcome), testSet$AlphScore), print.auc = TRUE, 
       col = "black", print.auc.y = .35)
  plot(roc(testSet$outcome, testSet$DEOGEN2_score_med), print.auc = TRUE, 
       col = "blue", print.auc.y = .25, add = TRUE)
  plot(roc(testSet$outcome, testSet$glm_AlphDeogen), print.auc = TRUE, 
       col = "red", print.auc.y = .15, add = TRUE)
  plot(roc(testSet$outcome, testSet$Alph_null), print.auc = TRUE, 
       col = "grey", print.auc.y = .05, add = TRUE)
  legend(0.2, 0.3, legend=c("AlphScore", "DEOGEN2", "AlphDeogen", "Alph_null"),
         col=c("black", "blue","red", "grey"), lty=1, cex=0.8)
  title(main = paste("Hold-out-set number", i))
  
  score_performance_tbl <- rbind(score_performance_tbl, 
                       tibble(Alph=roc(testSet$outcome, testSet$AlphScore)$auc, 
                              Alph_null=roc(testSet$outcome, testSet$Alph_null)$auc, 
                              CADD=roc(testSet$outcome, testSet$CADD_raw)$auc, 
                              REVEL=roc(testSet$outcome, testSet$REVEL_score)$auc,
                              DEOGEN2=roc(testSet$outcome, testSet$DEOGEN2_score_med)$auc, 
                              AlphCadd=roc(testSet$outcome, testSet$glm_AlphCadd)$auc, 
                              AlphDeogen=roc(testSet$outcome, testSet$glm_AlphDeogen)$auc, 
                              AlphRevel=roc(testSet$outcome, testSet$glm_AlphRevel)$auc, 
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

meanSemAUCs<-rbind(score_mean,
                   score_sd,
                   score_se)

meanSemAUCs$meaning<-c(
  "mean",
  "sd",
  "sem"
)

write_tsv(x=meanSemAUCs, 
          file="meanSemAUCs.tsv")

score_performance_tbl_spread<-score_performance_tbl%>%
  select(-num_train, -num_test)%>%
  gather(key="method")%>%
  arrange(value)

score_performance_tbl_spread<-score_performance_tbl_spread%>%
  mutate(method=factor(method, levels=c("Alph_null","Alph","CADD","AlphCadd","DEOGEN2","AlphDeogen","REVEL","AlphRevel")))

plot_aucs_ClinVar<-ggplot(score_performance_tbl_spread, aes(x=method, y=value))+
  stat_summary(fun.y = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.4)+
  geom_jitter(width=0.15, color="blue4")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=10),
        axis.text.y = element_text(color="black", size=10))+
  coord_cartesian(ylim=c(0.5,1)) + 
  labs(x = "")+
  labs(y = "AUC (ClinVar test set)", size=12)

plot_aucs_ClinVar
ggsave(filename= "plot_aucs_ClinVar.pdf", plot=plot_aucs_ClinVar, height=5, width=4)

pairwTTest<-pairwise.t.test(score_performance_tbl_spread$value, 
                     score_performance_tbl_spread$method, 
                     paired=TRUE, 
                     p.adj = "bonf")


