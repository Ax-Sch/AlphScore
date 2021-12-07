library(tidyverse)
library(pROC)

test_dataset2<-read_csv("data/prediction/base_model_test_dataset2.csv.gz")
validation_dataset<-read_csv("data/prediction/base_model_validation_set.csv.gz")

test_dataset2$multi<-!grepl(fixed("F1-"), test_dataset2$pdb_file)

ggplot(test_dataset2)+
  geom_point(aes(x=predicted, y=CADD_raw, color=multi), alpha=0.3)+
  facet_wrap(~outcome)


ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted, fill=outcome==1))+
  facet_wrap(~multi, scales = "free_y")

ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted, fill=outcome==1))+
  facet_wrap(~protein_mean_b_factor>80)

ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted, fill=outcome==1))+
  facet_wrap(~protein_length>1000 )

ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted, fill=outcome==1))+
  facet_wrap(~pLI>0.5 )


ggplot(test_dataset2)+
  geom_histogram(aes(x=CADD_raw, fill=outcome==1))


validation_dataset<-validation_dataset[is.na(validation_dataset$gnomAD_genomes_AC) & is.na(validation_dataset$gnomAD_genomes_AC),]

ggplot(validation_dataset)+
  geom_point(aes(x=CADD_raw, y=Score_exp), alpha=0.3)

ggplot(validation_dataset)+
  geom_point(aes(x=predicted, y=Score_exp), alpha=0.3)






cor.test(validation_dataset$Score_exp, validation_dataset$predicted)

validation_dataset$outcome<-validation_dataset$Score_exp<0.5
ggplot(validation_dataset)+
  geom_histogram(aes(x=predicted2, fill=outcome==1))

validation_dataset$outcome<-validation_dataset$Score_exp<0.5
ggplot(validation_dataset)+
  geom_histogram(aes(x=CADD_raw, fill=outcome==1))

validation_dataset$outcome<-validation_dataset$Score_exp<0.3

ggplot(validation_dataset)+
  geom_histogram(aes(x=predicted, fill=outcome==1))

for_comp<-validation_dataset[,c("Score_exp","predicted")]
for_comp<-for_comp %>% mutate(Score_exp_05=Score_exp<0.3, predicted_05=predicted>0.5)
table(for_comp[,c("Score_exp_05","predicted_05")])

roc_rose <- plot(roc(validation_dataset$outcome, validation_dataset$predicted), print.auc = TRUE, col = "red")

roc_rose <- plot(roc(validation_dataset$outcome, validation_dataset$predicted2), print.auc = TRUE, col = "red")

roc_rose <- plot(roc(validation_dataset$outcome, validation_dataset$CADD_raw), print.auc = TRUE, col = "red")

roc_rose <- plot(roc(validation_dataset$outcome, validation_dataset$REVEL_score), print.auc = TRUE, col = "red")

validation_dataset$REVEL_plus_ALPH<-validation_dataset$predicted_Alph + validation_dataset$REVEL_score

#str_split(validation_dataset$DEOGEN2_score,";")


spearmans_joined<-tibble()
for (un_ID in unique(validation_dataset$Uniprot_acc_split)){
  temp_values_joined<-validation_dataset %>% filter(Uniprot_acc_split==un_ID)
  for (dms in unique(temp_values_joined$DMS)){
    temp_inner_loop_values_joined<-temp_values_joined %>% filter(DMS==dms)
    for (score in c("predicted_Alph", "predicted_glm", "CADD_raw", "REVEL_score", "REVEL_plus_ALPH")){
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

spearmans_joined %>% 
  group_by(method) %>%
  summarise(mean(abs(spearm)))

spearmans_joined_spread<-spearmans_joined%>% spread(method, spearm)



for (gene_s in unique(validation_dataset$gene)){
  gene_filtered_set<-validation_dataset%>%
    filter(gene==gene_s)
  cadd_est<-cor.test(gene_filtered_set$Score_exp, gene_filtered_set$CADD_raw, method = "spearman")$estimate
  predicted_est<-cor.test(gene_filtered_set$Score_exp, gene_filtered_set$predicted, method = "spearman")$estimate
  predicted2_est<-cor.test(gene_filtered_set$Score_exp, gene_filtered_set$predicted2, method = "spearman")$estimate
  revel_est<-cor.test(gene_filtered_set$Score_exp, gene_filtered_set$REVEL_score, method = "spearman")$estimate
  revel_predicted_est<-cor.test(gene_filtered_set$Score_exp, gene_filtered_set$REVEL_score + gene_filtered_set$predicted, method = "spearman")$estimate
  
  print(paste(cadd_est, predicted_est, predicted2_est, revel_est, revel_predicted_est))
}
