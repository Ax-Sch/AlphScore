library(tidyverse)
library(pROC)

test_dataset2<-read_csv("data/prediction/base_model_test_dataset2.csv.gz")
validation_dataset<-read_csv("data/prediction/base_model_validation_set.csv.gz")
validation_dataset2<-read_csv("data/validation_set/values_joined_prepro_AlphScore.csv.gz")

validation_dataset<-validation_dataset2
test_dataset2$multi<-!grepl(fixed("F1-"), test_dataset2$pdb_file)

ggplot(test_dataset2)+
  geom_point(aes(x=predicted_Alph, y=CADD_raw, color=multi), alpha=0.3)+
  facet_wrap(~outcome)


ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted_Alph, fill=outcome==1))+
  facet_wrap(~multi, scales = "free_y")

ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted_Alph, fill=outcome==1))+
  facet_wrap(~protein_mean_b_factor>80)

ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted_Alph, fill=outcome==1))+
  facet_wrap(~protein_length>1000 )

ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted_Alph, fill=outcome==1))+
  facet_wrap(~pLI>0.5 )


ggplot(test_dataset2)+
  geom_histogram(aes(x=CADD_raw, fill=outcome==1))


validation_dataset<-validation_dataset[is.na(validation_dataset$gnomAD_genomes_AC) & is.na(validation_dataset$gnomAD_genomes_AC),]

ggplot(validation_dataset)+
  geom_point(aes(x=CADD_raw, y=DMS_val), alpha=0.3)

ggplot(validation_dataset)+
  geom_point(aes(x=predicted_Alph, y=DMS_val), alpha=0.3)

cor.test(validation_dataset$DMS_val, validation_dataset$AlphScore)


for_comp<-validation_dataset[,c("DMS_val","AlphScore")]
for_comp<-for_comp %>% mutate(DMS_val_05=DMS_val<0.3, predicted_05=AlphScore>0.5)
table(for_comp[,c("DMS_val_05","predicted_05")])

roc_rose <- plot(roc(validation_dataset$outcome, validation_dataset$AlphScore), print.auc = TRUE, col = "red")

roc_rose <- plot(roc(validation_dataset$outcome,  validation_dataset$CADD_raw+10*validation_dataset$AlphScore), print.auc = TRUE, col = "red")

roc_rose <- plot(roc(validation_dataset$outcome, validation_dataset$CADD_raw), print.auc = TRUE, col = "red")

roc_rose <- plot(roc(validation_dataset$outcome, validation_dataset$REVEL_score), print.auc = TRUE, col = "red")

validation_dataset$REVEL_plus_ALPH<-validation_dataset$AlphScore + validation_dataset$REVEL_score
validation_dataset$CADD_plus_ALPH<-validation_dataset$AlphScore*10 + validation_dataset$CADD_raw

validation_dataset$DEOGEN2_score<-str_replace(validation_dataset$DEOGEN2_score, fixed(";."), "")
validation_dataset$DEOGEN2_score<-(str_replace(validation_dataset$DEOGEN2_score, fixed(".;"), ""))

abc<-as.tibble(str_split(validation_dataset$DEOGEN2_score,";", simplify = TRUE))

abc %>% as.data.frame %>% rowwise() 

spearmans_joined<-tibble()
for (un_ID in unique(validation_dataset$Uniprot_acc_split)){
  temp_values_joined<-validation_dataset %>% filter(Uniprot_acc_split==un_ID)
  for (dms in unique(temp_values_joined$DMS)){
    temp_inner_loop_values_joined<-temp_values_joined %>% filter(DMS==dms)
    for (score in c("AlphScore", "CADD_raw", "CADD_plus_ALPH", "REVEL_score", "REVEL_plus_ALPH")){
    spearmans<-tibble(
      gene=unique(temp_inner_loop_values_joined$gene_dms), 
      UP_ID=un_ID, 
      DMS=dms, 
      spearm=cor.test(temp_inner_loop_values_joined$DMS_val, 
                        as.vector(unlist(temp_inner_loop_values_joined %>% select(one_of(score)))) , method = "spearman")$estimate,
    method=score)
        spearmans_joined=rbind(spearmans_joined,spearmans)
      }
  }
}

spearmans_joined_spread<-spearmans_joined%>% spread(method, spearm)

spearmans_joined %>% 
  group_by(method) %>%
  summarise(mean(abs(spearm)))
