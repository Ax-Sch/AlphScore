library(tidyverse)
library(pROC)
library(gridExtra)
library(optparse)
set.seed(1)

option_list = list(
  make_option(c("-t", "--variants"), type="character", default="data/prediction/base_model_variants.csv.gz", 
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
get_DEOGEN2_score<-function(DEOGEN2_score) {
  DEOGEN2_score<-str_replace(DEOGEN2_score, fixed(";."), "")
  DEOGEN2_score<-(str_replace(DEOGEN2_score, fixed(".;"), ""))
  median_score<-function(score_field){
    score_as_list<-str_split(score_field, ";", simplify = TRUE)
    score_as_list[score_as_list=="."]<-NA
    score_as_list<-as.double(score_as_list)
    return(median(score_as_list, na.rm=TRUE) )
  }
  DEOGEN2_score_med<-sapply(1:length(DEOGEN2_score), function(x) { median_score(DEOGEN2_score[x])})
  return(DEOGEN2_score_med)
}
variants$DEOGEN2_score_med<-get_DEOGEN2_score(variants$DEOGEN2_score)
variants$AlphScore<-variants$predicted_Alph

### combine scores on interim dataset:
interim_dataset<-variants %>% 
  filter(filter(CVinterim_no21_18_no_gnomad==TRUE,
                gnomadSet==FALSE))

model_glm_AC <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_dataset %>% dplyr::select(outcome, AlphScore, CADD_raw) %>%
                      filter(complete.cases(.)))

model_glm_AR <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_dataset %>% dplyr::select(outcome, AlphScore, REVEL_score) %>%
                      filter(complete.cases(.)))

model_glm_RC <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_dataset %>% dplyr::select(outcome, CADD_raw, REVEL_score) %>%
                      filter(complete.cases(.)))

model_glm_ARC <- glm(outcome ~ . , family=binomial(link='logit'),
                     data=interim_dataset %>% dplyr::select(outcome, AlphScore, CADD_raw, REVEL_score) %>%
                       filter(complete.cases(.)))

model_glm_AD <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_dataset %>% dplyr::select(outcome, AlphScore, DEOGEN2_score_med) %>%
                      filter(complete.cases(.)))

model_glm_CD <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_dataset %>% dplyr::select(outcome, CADD_raw, DEOGEN2_score_med) %>%
                      filter(complete.cases(.)))

model_glm_DR <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_dataset %>% dplyr::select(outcome, REVEL_score, DEOGEN2_score_med) %>%
                      filter(complete.cases(.)))

model_glm_ADR <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_dataset %>% dplyr::select(outcome, AlphScore, REVEL_score, DEOGEN2_score_med) %>%
                      filter(complete.cases(.)))

model_glm_ACD <- glm(outcome ~ . , family=binomial(link='logit'),
                     data=interim_dataset %>% dplyr::select(outcome, AlphScore, CADD_raw, DEOGEN2_score_med) %>%
                       filter(complete.cases(.)))

model_glm_CDR <- glm(outcome ~ . , family=binomial(link='logit'),
                     data=interim_dataset %>% dplyr::select(outcome, CADD_raw, DEOGEN2_score_med, REVEL_score) %>%
                       filter(complete.cases(.)))

validation_dataset$DEOGEN2_score_med<-get_DEOGEN2_score(validation_dataset$DEOGEN2_score)


validation_dataset$predicted_glm_AC<-predict(model_glm_AC, validation_dataset)
validation_dataset$predicted_glm_AR<-predict(model_glm_AR, validation_dataset)
validation_dataset$predicted_glm_RC<-predict(model_glm_RC, validation_dataset)
validation_dataset$predicted_glm_ARC<-predict(model_glm_ARC, validation_dataset)
validation_dataset$predicted_glm_AD<-predict(model_glm_AD, validation_dataset)
validation_dataset$predicted_glm_CD<-predict(model_glm_CD, validation_dataset)
validation_dataset$predicted_glm_DR<-predict(model_glm_DR, validation_dataset)
validation_dataset$predicted_glm_ADR<-predict(model_glm_ADR, validation_dataset)
validation_dataset$predicted_glm_ACD<-predict(model_glm_ACD, validation_dataset)
validation_dataset$predicted_glm_CDR<-predict(model_glm_CDR, validation_dataset)

###### validation Dataset


validation_dataset<-validation_dataset %>%
  filter(is.na(gnomAD_genomes_AC) & is.na(gnomAD_genomes_AC))

check_corr_CADD<-ggplot(validation_dataset)+
  geom_point(aes(x=CADD_raw, y=DMS_val), alpha=0.1)+
  facet_wrap(~ paste(DMS, gene_dms), scales="free")
print(check_corr_CADD)
ggsave(filename="check_corr_CADD.pdf", plot=check_corr_CADD)

check_corr_Alph<-ggplot(validation_dataset)+
  geom_point(aes(x=AlphScore, y=DMS_val), alpha=0.1)+
  facet_wrap(~ paste(DMS, gene_dms), scales="free")
print(check_corr_Alph)
ggsave(filename="check_corr_Alph.pdf", plot=check_corr_Alph)



SCORES<-c("AlphScore", "CADD_raw", "predicted_glm_AC", 
          "REVEL_score", "predicted_glm_AR", "predicted_glm_RC", 
          "DEOGEN2_score_med", "predicted_glm_AD", "predicted_glm_CD", 
          "predicted_glm_ACD", "predicted_glm_DR","predicted_glm_ADR",
          "predicted_glm_CDR")
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


friedman.test(abs(spearmans_joined$spearm), 
              spearmans_joined$method, 
              paste(spearmans_joined$DMS, spearmans_joined$UP_ID))

pairwise.wilcox.test(abs(spearmans_joined$spearm), 
                     spearmans_joined$method, 
                     paired=TRUE, 
                     p.adj = "bonf")
