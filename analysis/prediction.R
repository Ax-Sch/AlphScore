###### 19.10.2021 todo:
### prefilter variants is okay
### build the best score possible with just using alphafold-derived features
### (interpretability)
### then possibly build a score with conservation etc. included
### check the distribution of variant exchanges in training vs. test set
### build delta scores for to_AS variables
### ? include a score that captures the AS-exchange?

# 04.11.2021 Todo
### CpG Varianten in evolutionär jüngeren Varianten, Deaminierung von phosphorylierten 
### verhindern, dass er die ausgehende Aminosäure als Feature entdeckt
### Cytoscape einmal ansehen
### stratified split für hold out.
### Tiefe der Bäume
### ggf. systematisch durchprobieren
### Ziel: hold out perfomane darf fallen; test performance steigt
### hold out set
### Hydrophobicity, Polarität
### Feature: wie homogen ist die Nachbarschaft
### was passiert mit der Feature importance, wenn man die Daten vorfiltert?
### Ziel: viele Feature tragen bei, extremer Abfall ist ein schlechtes Zeichen

### 23.11. : 
### Mutationen matchen zwischen den zwei Labels 
### Deep Mutational Scan testen
### 2018-2021 evtl. kritisch; evtl zufällig besser
### Snakemake workflow teilbar machen
### Daten vorberechnen
### Ggf. EVE ansehen, wegen Scan-Daten
### Martin Kircher schickt Scan-Daten; MAVE-Datenbank
### Genomic Medicine 



library(tidyverse)
library(pROC)
library("readxl")
library(caret)

args=c("variants_preprocessed_recalibrated_v2.csv.gz", "available_colnames_W_surr.xlsx", 
       "0_0_xgboost_300_6_0_0.05_0.8_0_0.999999_90_5_0", 
       "0", "0", "extra", "2000", "0", "0", "0", "0", "0", "0.999999", "90", "0", "0")

args=c("variants_preprocessed_recalibrated_v2.csv.gz", "available_colnames_W_surr.xlsx", 
       "0_0_xgboost_300_6_0_0.05_0.8_0_0.999999_90_5_0", 
       "0", "0", "xgboost", "300", "6", "0", "0.05", "0.8", "0", "0.999999", "90", "0", "0")

args=c("variants_preprocessed_recalibrated_v2_strict.csv.gz", "available_colnames_W_surr.xlsx", 
       "0_0_xgboost_300_6_0_0.05_0.8_0_0.999999_90_5_0", 
       "0", "0", "xgboost", "300", "6", "0", "0.05", "0.8", "0", "0.999999", "90", "0", "0")

args=c("variants_preprocessed_v2.csv.gz", "available_colnames_W_surr.xlsx", 
       "0_0_xgboost_300_6_0_0.05_0.8_0_0.999999_90_5_0", 
       "0", "0", "xgboost", "300", "6", "0", "0.05", "0.8", "0", "0.999999", "90", "0", "0")



args=c("variants_preprocessed_recalibrated_v2.csv.gz", "available_colnames_W_surr.xlsx", 
       "0_0_xgboost_300_6_0_0.05_0.8_0_0.999999_90_5_0", 
       "0", "0", "randomforest", "2000", "5", "0", "0", "0", "0", "0.999999", "90", "0", "0")


args=c("variants_preprocessed_recalibrated_v2_strict.csv.gz", "available_colnames_W_surr.xlsx", 
       "0_0_xgboost_300_6_0_0.05_0.8_0_0.999999_90_5_0", 
       "0", "0", "randomforest", "2000", "5", "0", "0", "0", "0", "0.999999", "90", "0", "0")

# nicht ausführen!:
args = commandArgs(trailingOnly=TRUE)

print(args)

csv_location=args[1]
excel_locataion=args[2]
prefix=args[3]
filter_factor_specific=as.double(args[4])
filter_factor_global=as.double(args[5])
method_pred=as.character(args[6])
num_trees_param=as.integer(args[7])
max_depth_param=as.integer(args[8])
lambda_param=as.double(args[9])
eta_param=as.double(args[10])
subsample_param=as.double(args[11])
min_node_param=as.integer(args[12])
cor_param=as.double(args[13])
b_factor_param=as.integer(args[14])
min_child_weight_param=as.integer(args[15])
gamma_param=as.integer(args[16])

pdf(file=paste0(prefix, "_RPlots.pdf"))
#pdf(NULL)


variants_org<-read_csv(csv_location, na=c(".","NA", NA))
write_tsv(as_tibble(colnames(variants_org)), file="available_colnames.tsv")
colnames_usage <- read_excel(excel_locataion, col_types="text")

#### experimental modifications of variables:
mean_pli<-mean(variants_org$pLI, na.rm = TRUE)
variants<-variants_org %>%
  mutate(pLI=ifelse(is.na(pLI), mean_pli,pLI))

sel_vars_to<-(colnames_usage %>% filter(!is.na(add_to_AS)))$value

toAS_properties<-variants  %>%
  filter(gnomadSet == 1, b_factor>b_factor_param)%>%
  group_by(from_AS) %>%
  dplyr::select(sel_vars_to)%>%
  summarize_all(mean)
colnames(toAS_properties)<-paste0(colnames(toAS_properties), "_toAS")


variants<-variants%>%
  left_join(toAS_properties, by=c("from_AS"="from_AS_toAS"))%>%
  mutate(CADD_raw_scaled=scale(CADD_raw))

variants[, colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])]<-
  variants[, sel_vars_to] - 
  variants[, colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])]


colnames_prediction<-c(colnames(toAS_properties %>% dplyr::select(-from_AS_toAS)), 
                       (colnames_usage %>% filter(!is.na(for_prediction)))$value)



train_dataset<-variants %>% 
    filter(!(pure_cv18_to_21_gene==TRUE) & gnomadSet==TRUE)%>%
    mutate(mean_revel_global=mean(REVEL_score, na.rm=TRUE)) %>%
    group_by(from_AS, to_AS, outcome)%>%
    mutate(mean_revel=mean(REVEL_score, na.rm=TRUE))%>%
    ungroup()%>%
    filter((REVEL_score>(mean_revel*filter_factor_specific)) | outcome==0)%>%
    filter((REVEL_score>(mean_revel_global*filter_factor_global)) | outcome==0)%>%
    dplyr::select(c(colnames_prediction, "outcome"))


#https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
zv <- apply(train_dataset, 2, function(x) length(unique(x)) == 1)
dfr <- train_dataset[, !zv]
n=length(colnames(dfr))
correlationMatrix <- cor(dfr,use="complete.obs")


cols_removed<-dfr[, -(findCorrelation(correlationMatrix, cutoff=cor_param))]
colnames_new<-colnames(cols_removed)
#colnames_new<-colnames_new[!colnames_new %in% remove_vars]
write_tsv(as_tibble(colnames_new), file="variables_used.tsv")

train_dataset<-train_dataset[,colnames_new] 
train_dataset<-train_dataset[complete.cases(train_dataset),]

set.seed(1)

test_dataset<-variants %>% 
  filter(clinvar_no_cv21to18_no_gnomad==TRUE)
test_dataset<-test_dataset[complete.cases(test_dataset[,c(colnames_prediction, "outcome", "CADD_raw")]), ]


ranger_fit<-function(xtrain_dataset){
  library(ranger)
  model1<-ranger(outcome ~ ., 
                 data=xtrain_dataset, 
                 importance="impurity", max.depth=max_depth_param,
                 num.trees = num_trees_param,
                 min.node.size = min_node_param
  ) 
  return(model1)
}

ranger_predict<-function(dataset, modelx){
  return(predict(modelx, dataset)$predictions)
}


extraT_fit<-function(xtrain_dataset){
  library(extraTrees)
  model1<-extraTrees(
                 x = as.matrix(xtrain_dataset %>% dplyr::select(-outcome)),
                 y =xtrain_dataset$outcome,
                 ntree = num_trees_param,
                 numThreads=7
  ) 
  return(model1)
}

extraT_predict<-function(dataset, modelx){
  
  return(
  predict(modelx, as.matrix(dataset %>% dplyr::select(colnames_new) %>% dplyr::select(-outcome) )) 
)
}



xgboost_fit<-function(xtrain_dataset){
  library(xgboost)
  n_ds<-nrow(xtrain_dataset)
  train_i<-sample(1:n_ds,as.integer(n_ds*0.8))
  train_ds<-xgb.DMatrix(data=as.matrix(xtrain_dataset[train_i,] %>% dplyr::select(-outcome)), label=xtrain_dataset$outcome[train_i] )
  #test_ds<-xgb.DMatrix(data=as.matrix(test_dataset %>% dplyr::select(colnames_new, -outcome)), label=test_dataset$outcome )
  
  test_ds<-xgb.DMatrix(data=as.matrix(xtrain_dataset[-train_i,] %>% dplyr::select(-outcome)), label=xtrain_dataset$outcome[-train_i] )
  watchlist=list(train=train_ds, test=test_ds)
  
  params <- list(max_depth = max_depth_param, subsample=subsample_param, 
                 eta=eta_param, lambda=lambda_param, alpha=0, 
                 gamma=gamma_param, min_child_weight=min_child_weight_param)
  model1<-xgb.train(
    data=train_ds,
    watchlist = watchlist,
    params=params,
    nrounds = num_trees_param,
    early_stopping_rounds = 10,
    #eval_metric="auc",
    objective = "binary:logistic",  # for regression models
    verbose = 1)
  
  return(model1)
}

xgboost_predict<-function(dataset, modelx){
  data_x<-xgb.DMatrix(data=as.matrix(dataset %>% dplyr::select(-outcome)))
  return(
    predict(modelx, data_x) 
  )
}

if (method_pred=="xgboost"){
fit_model<-xgboost_fit
predict_model<-xgboost_predict
} else if (method_pred=="randomforest") {
  fit_model<-ranger_fit
  predict_model<-ranger_predict
}else {
  fit_model<-extraT_fit
  predict_model<-extraT_predict
}


model1<-fit_model(train_dataset %>% dplyr::select(colnames_new))

train_dataset$predicted<-predict_model(train_dataset %>% dplyr::select(colnames_new), model1)
roc_rose <- plot(roc(train_dataset$outcome, train_dataset$predicted), print.auc = TRUE, col = "red")
auc_train=roc_rose$auc


test_dataset$predicted<-predict_model(test_dataset %>% dplyr::select(colnames_new), model1)

roc_rose <- plot(roc(test_dataset$outcome, test_dataset$CADD_raw), print.auc = TRUE, col = "red")
roc_rose <- plot(roc(test_dataset$outcome, test_dataset$predicted), print.auc = TRUE, 
                 col = "green", print.auc.y = .4, add = TRUE)


auc1=roc_rose$auc



model2 <- glm(outcome ~ . , family=binomial(link='logit'),
              data=test_dataset %>% dplyr::select(outcome, predicted, CADD_raw) %>%
                filter(complete.cases(.)))

#library(e1071)
#test_dataset$predicted_scaled<-scale(test_dataset$predicted)

#classifierR = svm(formula = outcome ~ predicted_scaled + CADD_raw_scaled,
#                 data = test_dataset,
#                kernel = 'linear')



test_dataset2<-variants %>% 
  filter(cv18_to_21_CV_test==TRUE)

test_dataset2<-test_dataset2[complete.cases(test_dataset2[,c(colnames_new, "outcome")]), ]
test_dataset2$predicted<-predict_model(test_dataset2 %>% dplyr::select(colnames_new), model1)
test_dataset2$predicted_scaled<-scale(test_dataset2$predicted)

test_dataset2$predicted2<-predict(model2, test_dataset2)
#test_dataset2$predicted3<-predict(classifierR, test_dataset2[,c("predicted_scaled", "CADD_raw_scaled")])



roc_rose <- plot(roc(test_dataset2$outcome, test_dataset2$CADD_raw), print.auc = TRUE, col = "red")
roc_rose <- plot(roc(test_dataset2$outcome, test_dataset2$predicted2), print.auc = TRUE, 
                 col = "blue", print.auc.y = .2, add = TRUE)
roc_rose <- plot(roc(test_dataset2$outcome, test_dataset2$predicted), print.auc = TRUE, 
                 col = "green", print.auc.y = .4, add = TRUE)
#roc_rose <- plot(roc(test_dataset2$outcome, test_dataset2$predicted3), print.auc = TRUE, 
#                 col = "black", print.auc.y = .6, add = TRUE)

print(roc_rose)

auc2=roc_rose$auc
save_tibble<-data.frame(auc1=auc1, auc2=auc2, condition=prefix, params=I(list(args)), auc_train=auc_train)

write.csv2(x=save_tibble, file=paste0(prefix, "_results.tsv"))


if (method_pred!="xgboost"){
##### CHECK properties of models
var_imp<-tibble(importance=as.vector(model1$variable.importance), variable=names(model1$variable.importance))

var_importance<-ggplot(var_imp, aes(x=reorder(variable,importance), y=importance,fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  ggtitle("Information Value Summary")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")

remove_vars<-(var_imp %>% arrange(importance))[1:40,]$variable

ggsave(filename=paste0(prefix,"_importance.pdf"), plot=var_importance, height=49)
}
ggplot(test_dataset2)+
  geom_point(aes(x=predicted, y=CADD_raw), alpha=0.3)+
  facet_wrap(~outcome)

