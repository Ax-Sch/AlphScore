library(tidyverse)
library(pROC)
library(readxl)
library(caret)
library(optparse)
set.seed(1)

option_list = list(
  make_option(c("-c", "--csv_location"), type="character", default="/media/axel/Dateien/Arbeit_Gen/alphafold2/data_from_xcat_v2/variants_preprocessed_recalibrated_v2.csv.gz", 
              help="csv.gz file"),
  make_option(c("-e", "--excel_locataion"), type="character", default="config/available_colnames_W_surr.xlsx", 
              help="Excel file listing columns to use"),
  make_option(c("-p", "--prefix"), type="character", default="abc", 
              help="Prefix for output"),
  make_option(c("-m", "--method_pred"), type="character", default="randomforest", 
              help="Prediction method [randomforest, xgboost, extratree]"),
  make_option(c("-n", "--num_trees_param"), type="integer", default=2000, 
              help="Number of trees to use / or number of rounds"),
  make_option(c("-x", "--max_depth_param"), type="integer", default=6, 
              help="Maximal depth of trees"),
  make_option(c("-l", "--lambda_param"), type="double", default=0, 
              help="Lambda value for xgboost"),
  make_option(c("-t", "--eta_param"), type="double", default=0.05, 
              help="Eta - learning rate, xgboost"),
  make_option(c("-s", "--subsample_param"), type="double", default=0.8, 
              help="Subsample parameter, xgboost"),
  make_option(c("-i", "--min_node_param"), type="integer", default=50, 
              help="Min Node parameter, xgboost"),
  make_option(c("-r", "--cor_param"), type="double", default=0.999999, 
              help="Just keep columns that have a lower correlation than this value, general"),
  make_option(c("-b", "--b_factor_param"), type="integer", default=90, 
              help="pLDDT value to filter positions that will be used to generate the average values of amino acids, general"),
  make_option(c("-d", "--min_child_weight_param"), type="integer", default=1, 
              help="min child weight, xgboost"),
  make_option(c("-g", "--gamma_param"), type="integer", default=90, 
              help="gamma value, xgboost"),
  make_option(c("-o", "--out_folder"), type="character", default="prediction", 
              help="gamma value, xgboost")
  
)

opt = parse_args(OptionParser(option_list=option_list))
#DEBUG:
#opt$csv_location="/media/axel/Dateien/Arbeit_Gen/alphafold2/data_from_xcat_v2/variants_preprocessed_recalibrated_v2.csv.gz"



variants<-read_csv(opt$csv_location, na=c(".","NA", NA))
colnames_usage <- read_excel(opt$excel_locataion, col_types="text")

work_dir=paste0("data/", opt$out_folder)
dir.create(work_dir)
setwd(work_dir)
pdf(file=paste0(opt$prefix, "_RPlots.pdf"))

sel_vars_to<-(colnames_usage %>% filter(!is.na(add_to_AS)))$value
toAS_properties<-variants  %>%
  filter(gnomadSet == 1, b_factor>opt$b_factor_param)%>%
  group_by(from_AS) %>%
  dplyr::select(all_of(sel_vars_to))%>%
  summarize_all(mean)
colnames(toAS_properties)<-paste0(colnames(toAS_properties), "_toAS")

variants<-variants%>%
  left_join(toAS_properties, by=c("from_AS"="from_AS_toAS"))

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
    dplyr::select(c(colnames_prediction, "outcome"))

#https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
zv <- apply(train_dataset, 2, function(x) length(unique(x)) == 1)
dfr <- train_dataset[, !zv]
n=length(colnames(dfr))
correlationMatrix <- cor(dfr,use="complete.obs")

cols_removed<-dfr[, -(findCorrelation(correlationMatrix, cutoff=opt$cor_param))]
colnames_new<-colnames(cols_removed)
write_tsv(as_tibble(colnames_new), file="variables_used.tsv")

train_dataset<-train_dataset[,colnames_new] 
train_dataset<-train_dataset[complete.cases(train_dataset),]

test_dataset<-variants %>% 
  filter(clinvar_no_cv21to18_no_gnomad==TRUE)
test_dataset<-test_dataset[complete.cases(test_dataset[,c(colnames_prediction, "outcome", "CADD_raw")]), ]


ranger_fit<-function(xtrain_dataset){
  library(ranger)
  model1<-ranger(outcome ~ ., 
                 data=xtrain_dataset, 
                 importance="impurity", max.depth=opt$max_depth_param,
                 num.trees = opt$num_trees_param,
                 min.node.size = opt$min_node_param
  ) 
  return(model1)
}
ranger_predict<-function(dataset, modelx){
  return(predict(modelx, dataset)$predictions) }

extraT_fit<-function(xtrain_dataset){
  library(extraTrees)
  model1<-extraTrees(
                 x = as.matrix(xtrain_dataset %>% dplyr::select(-outcome)),
                 y =xtrain_dataset$outcome,
                 ntree = opt$num_trees_param,
                 numThreads=7
  ) 
  return(model1)
}
extraT_predict<-function(dataset, modelx){
  return(
  predict(modelx, as.matrix(dataset %>% dplyr::select(all_of(colnames_new)) %>% dplyr::select(-outcome)))   )}

xgboost_fit<-function(xtrain_dataset){
  library(xgboost)
  n_ds<-nrow(xtrain_dataset)
  train_i<-sample(1:n_ds,as.integer(n_ds*0.8))
  train_ds<-xgb.DMatrix(data=as.matrix(xtrain_dataset[train_i,] %>% dplyr::select(-outcome)), label=xtrain_dataset$outcome[train_i] )
  #test_ds<-xgb.DMatrix(data=as.matrix(test_dataset %>% dplyr::select(colnames_new, -outcome)), label=test_dataset$outcome )
  test_ds<-xgb.DMatrix(data=as.matrix(xtrain_dataset[-train_i,] %>% dplyr::select(-outcome)), label=xtrain_dataset$outcome[-train_i] )
  watchlist=list(train=train_ds, test=test_ds)
  
  params <- list(max_depth = opt$max_depth_param, subsample=opt$subsample_param, 
                 eta=opt$eta_param, lambda=opt$lambda_param, alpha=0, 
                 gamma=opt$gamma_param, min_child_weight=opt$min_child_weight_param)
  model1<-xgb.train(
    data=train_ds,
    watchlist = watchlist,
    params=params,
    nrounds = opt$num_trees_param,
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

if (opt$method_pred=="xgboost"){
fit_model<-xgboost_fit
predict_model<-xgboost_predict
} else if (opt$method_pred=="randomforest") {
  fit_model<-ranger_fit
  predict_model<-ranger_predict
}else {
  fit_model<-extraT_fit
  predict_model<-extraT_predict
}



model1<-fit_model(train_dataset %>% dplyr::select(all_of(colnames_new)))

train_dataset$predicted<-predict_model(train_dataset %>% dplyr::select(all_of(colnames_new)), model1)
roc_rose <- plot(roc(train_dataset$outcome, train_dataset$predicted), print.auc = TRUE, col = "red")
auc_train=roc_rose$auc


test_dataset$predicted<-predict_model(test_dataset %>% dplyr::select(all_of(colnames_new)), model1)

roc_rose <- plot(roc(test_dataset$outcome, test_dataset$CADD_raw), print.auc = TRUE, col = "red")
roc_rose <- plot(roc(test_dataset$outcome, test_dataset$predicted), print.auc = TRUE, 
                 col = "green", print.auc.y = .4, add = TRUE)


auc1=roc_rose$auc



model2 <- glm(outcome ~ . , family=binomial(link='logit'),
              data=test_dataset %>% dplyr::select(outcome, predicted, CADD_raw) %>%
                filter(complete.cases(.)))


test_dataset2<-variants %>% 
  filter(cv18_to_21_CV_test==TRUE)

test_dataset2<-test_dataset2[complete.cases(test_dataset2[,c(colnames_new, "outcome")]), ]
test_dataset2$predicted<-predict_model(test_dataset2 %>% dplyr::select(all_of(colnames_new)), model1)
test_dataset2$predicted_scaled<-scale(test_dataset2$predicted)

test_dataset2$predicted2<-predict(model2, test_dataset2)


roc_rose <- plot(roc(test_dataset2$outcome, test_dataset2$CADD_raw), print.auc = TRUE, col = "red")
roc_rose <- plot(roc(test_dataset2$outcome, test_dataset2$predicted2), print.auc = TRUE, 
                 col = "blue", print.auc.y = .2, add = TRUE)
roc_rose <- plot(roc(test_dataset2$outcome, test_dataset2$predicted), print.auc = TRUE, 
                 col = "green", print.auc.y = .4, add = TRUE)
print(roc_rose)

auc2=roc_rose$auc
save_tibble<-data.frame(auc1=auc1, auc2=auc2, condition=opt$prefix, params=I(list(opt)), auc_train=auc_train)

write.csv2(x=save_tibble, file=paste0(opt$prefix, "_results.tsv"))


if (opt$method_pred=="randomforest"){
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

ggsave(filename=paste0(opt$prefix,"_importance.pdf"), plot=var_importance, height=49)
}
ggplot(test_dataset2)+
  geom_point(aes(x=predicted, y=CADD_raw), alpha=0.3)+
  facet_wrap(~outcome)

