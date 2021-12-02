library(tidyverse)
library(pROC)

test_dataset2<-read_csv("data/prediction/abc_test_dataset2.csv.gz")
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
