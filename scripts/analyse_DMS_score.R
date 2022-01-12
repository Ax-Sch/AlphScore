library(tidyverse)
library(pROC)
library(gridExtra)
library(optparse)
source("scripts/existing_scores_glm_functions.R")
set.seed(1)

option_list = list(
  make_option(c("-t", "--variants"), type="character", default="data/prediction_final/pre_final_model_regular_variants.csv.gz", 
              help="csv.gz file, test dataset"),
  make_option(c("-v", "--validation_set"), type="character", default="data/validation_set/validation_set_w_AlphScore.csv.gz", 
              help="Validation set"),
  make_option(c("-o", "--out_folder"), type="character", default="data/analyse_score", 
              help="Output folder"))

opt = parse_args(OptionParser(option_list=option_list))

variants<-read_csv(opt$variants)
validation_dataset<-read_csv(opt$validation_set)

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

# get DEOGEN2 Scores:
variants$pos_in_VEP_and_Uniprot<-get_index_col(variants)
variants$DEOGEN2_score_med<-unlist_score(variants$DEOGEN2_score, variants$pos_in_VEP_and_Uniprot)

validation_dataset$pos_in_VEP_and_Uniprot<-get_index_col(validation_dataset)
validation_dataset$DEOGEN2_score_med<-unlist_score(validation_dataset$DEOGEN2_score, validation_dataset$pos_in_VEP_and_Uniprot)

variants$AlphScore<-variants$predicted_Alph
### combine scores on interim dataset:
interim_dataset<-variants %>% 
  filter(clinvar_interim_test==TRUE)


set_of_models<-fit_set_of_models(interim_dataset)
validation_dataset<-predict_set_of_models(set_of_models, validation_dataset)

###### validation Dataset

variants<-variants%>% 
  mutate(ID=paste(`#chr`, `pos(1-based)`, ref, alt, sep=":"))

validation_dataset<-validation_dataset %>%
  mutate(ID=paste(`#chr`, `pos(1-based)`, ref, alt, sep=":"))%>%
  filter(!ID %in% variants$ID)

SCORES<-c("AlphScore", "CADD_raw", "glm_AlphCadd", 
          "REVEL_score", "glm_AlphRevel", "glm_RevelCadd", 
          "DEOGEN2_score_med", "glm_AlphDeogen", "glm_CaddDeogen", 
          "glm_AlphCaddDeogen", "glm_DeogenRevel","glm_AlphDeogenRevel",
          "glm_CaddDeogenRevel")
spearmans_joined<-tibble()
for (un_ID in unique(validation_dataset$Uniprot_acc_split)){
  temp_values_joined<-validation_dataset %>% 
    filter(Uniprot_acc_split==un_ID)
  for (dms in unique(temp_values_joined$DMS)){
    temp_inner_loop_values_joined<-temp_values_joined %>% 
      filter(DMS==dms)
    for (score in SCORES){
     score_DMS_tibble<-tibble(score_val=unlist(temp_inner_loop_values_joined %>% select(one_of(score))),
                              DMS_val=temp_inner_loop_values_joined$DMS_val)
     score_DMS_tibble<-score_DMS_tibble %>% filter(complete.cases(.))
      spearmans<-tibble(
      gene=unique(temp_inner_loop_values_joined$gene_dms), 
      UP_ID=un_ID, 
      DMS=dms, 
      spearm=ifelse(nrow(score_DMS_tibble)>2,
      cor.test(score_DMS_tibble$DMS_val, score_DMS_tibble$score_val,
               method = "spearman", na.action="na.ommit")$estimate,NA),
      method=score,
      number_of_vals=nrow(score_DMS_tibble))
        spearmans_joined=rbind(spearmans_joined,spearmans)
      }
  }
}

prot_info_for_join<-validation_dataset %>% 
  distinct(Uniprot_acc_split, protein_length, protein_mean_b_factor)

bad_UP_ids<-(spearmans_joined %>% 
               filter(number_of_vals<21))$UP_ID

spearmans_joined<-spearmans_joined %>% 
  mutate(method=factor(method, levels= SCORES))%>%
  filter(! UP_ID %in% bad_UP_ids)%>%
  left_join(prot_info_for_join, by=c("UP_ID"="Uniprot_acc_split"))


spearmans_joined_spread<-spearmans_joined%>% 
  spread(method, spearm)

spearmans_joined_summarized<-spearmans_joined %>% 
  group_by(method, protein_length, protein_mean_b_factor) %>%
  summarise(mean_abs_correlation=mean(abs(spearm), na.rm = TRUE))%>%
  mutate(method=factor(method, levels= SCORES)) 

plot_spearmans<-ggplot(spearmans_joined, aes(x=method, y=abs(spearm)))+
  stat_summary(fun.y = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.4)+
  geom_jitter(width=0.15)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1))

plot_spearmans
ggsave(filename= "spearman_plot.pdf", plot=plot_spearmans, height=8, width=6)

plot_spearmans<-ggplot(spearmans_joined, aes(x=method, y=abs(spearm)))+
  stat_summary(fun.y = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.4)+
  geom_jitter(width=0.15, aes(color=protein_mean_b_factor, fill=protein_length), size=2, shape=21)+
  scale_color_gradient(low="red", high="blue")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1))

plot_spearmans
ggsave(filename= "spearman_plot_colored.pdf", plot=plot_spearmans, height=8, width=8)

table_spearm<-spearmans_joined %>% 
  mutate(abs_spearman=abs(spearm))%>%
  group_by(method)%>%
  summarise(mean_spearm=mean(abs_spearman), sd=sd(abs_spearman), num=n(), sem=sd(abs_spearman)/sqrt(n())) %>%
  arrange(-mean_spearm)

print(table_spearm)

write_tsv(x = table_spearm, file = "table_spearm.tsv")

friedman.test(abs(spearmans_joined$spearm), 
              spearmans_joined$method, 
              paste(spearmans_joined$DMS, spearmans_joined$UP_ID))

pairwise.wilcox.test(abs(spearmans_joined$spearm), 
                     spearmans_joined$method, 
                     paired=TRUE, 
                     p.adj = "bonf")

