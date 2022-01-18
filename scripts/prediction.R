library(tidyverse)
library(pROC)
library(readxl)
library(caret)
library(optparse)
source("scripts/existing_scores_glm_functions.R")

set.seed(1)

option_list = list(
  make_option(c("-c", "--csv_location"), type="character", default="/media/axel/Dateien/Arbeit_Gen/alphafold2/git_version/AlphScore/data/train_testset1/subsampled_vars.csv", 
              help="csv.gz file"),
  make_option(c("-e", "--excel_location"), type="character", default="resources/available_colnames_regular.xlsx", 
              help="Excel file listing columns to use"),
  make_option(c("-p", "--prefix"), type="character", default="base_model", 
              help="Prefix for output"),
  make_option(c("-m", "--method_pred"), type="character", default="randomforest", 
              help="Prediction method [randomforest, xgboost, extratree]"),
  make_option(c("-n", "--num_trees_param"), type="integer", default=200, 
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
              help="randomforest: min.node.size, xgboost: min node parameter, extratree: nodesize"),
  make_option(c("-r", "--cor_param"), type="double", default=0.999999, 
              help="Just keep columns that have a lower correlation than this value, general"),
  make_option(c("-b", "--b_factor_param"), type="integer", default=20, 
              help="pLDDT value to filter positions that will be used to generate the average values of amino acids, general"),
  make_option(c("-d", "--min_child_weight_param"), type="integer", default=0, 
              help="min child weight, xgboost"),
  make_option(c("-g", "--gamma_param"), type="integer", default=0, 
              help="gamma value, xgboost"),
  make_option(c("-o", "--out_folder"), type="character", default="data/prediction", 
              help="name of folder to store output"),
  make_option(c("-w", "--write_dataset"), type="logical", default=FALSE, 
              help="Write predictions of test-dataset to file"),
  make_option(c("-f", "--full_model"), type="logical", default=FALSE, #data/preprocess/validation_set.csv.gzpreprocessed.csv.gz
              help="generate a full model"),
  make_option(c("-a", "--importance"), type="character", default="impurity", 
              help="variable importance mode of the ranger model"),
  make_option(c("-y", "--write_model"), type="logical", default=FALSE, #data/preprocess/validation_set.csv.gzpreprocessed.csv.gz
              help="write the model parameters to a file"),
  make_option(c("-k", "--k_fold_cross_val"), type="logical", default=FALSE,
              help="activate k-fold cross validation of the training and test data set")
  
)

opt = parse_args(OptionParser(option_list=option_list))

save_tibble<-data.frame()

variants<-read_csv(opt$csv_location, na=c(".","NA", NA))
colnames_usage <- read_excel(opt$excel_location, col_types="text")

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

pdf(file=paste0(opt$prefix, "_RPlots.pdf"))

# calculate average values of add_to_AS column for amino acids of training set
sel_vars_to<-(colnames_usage %>% filter(!is.na(add_to_AS)))$value
toAS_properties<-variants  %>%
  filter(gnomadSet == 1, b_factor>opt$b_factor_param)%>%
  dplyr::select(all_of(sel_vars_to), from_AS)%>%
  filter(complete.cases(.))%>%
  group_by(from_AS) %>%
  summarize_all(mean)
colnames(toAS_properties)<-paste0(colnames(toAS_properties), "_toAS")

# add these toAS values to the data set
variants<-variants%>%
  left_join(toAS_properties, by=c("from_AS"="from_AS_toAS"))
colnames_toAS<-colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])
variants[, colnames_toAS]<-  variants[, sel_vars_to] - 
  variants[, colnames_toAS]

# remove correlated columns in the data set
colnames_prediction<-c(colnames_toAS, 
                       (colnames_usage %>% 
                          filter(!is.na(for_prediction)))$value)

variants_pred<-variants %>% 
  dplyr::select(c(colnames_prediction, "outcome"))

#https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
zv <- apply(variants_pred, 2, function(x) length(unique(x)) == 1)
dfr <- variants_pred[, !zv]
n=length(colnames(dfr))
correlationMatrix <- cor(dfr,use="complete.obs")

set.seed(1)
# prune correltated variables 
df_colsCorrRemoved<-dfr[, -(findCorrelation(correlationMatrix, cutoff=opt$cor_param))]
colsCorrRemoved<-colnames(df_colsCorrRemoved)
rm(df_colsCorrRemoved)
rm(dfr)

variants<-variants[complete.cases(variants[,c(colsCorrRemoved, "outcome", "CADD_raw")]), ]
gc()


ranger_fit<-function(xtrain_dataset){
  library(ranger)
  model1<-ranger(outcome ~ ., 
                 data=xtrain_dataset, 
                 importance=opt$importance, 
                 max.depth=opt$max_depth_param,
                 num.trees = opt$num_trees_param,
                 min.node.size = opt$min_node_param
  ) 
  return(model1)
}
ranger_predict<-function(dataset, modelx){
  return(predict(modelx, dataset)$predictions) }

extraT_fit<-function(xtrain_dataset){
  library(ranger)
  model1<-ranger(outcome ~ ., 
                 data=xtrain_dataset, 
                 importance=opt$importance, 
                 splitrule="extratrees",
                 max.depth=opt$max_depth_param,
                 num.trees = opt$num_trees_param,
                 min.node.size = opt$min_node_param
  ) 
  return(model1)
}
ranger_save<-function(modelx, filename){
  saveRDS(modelx, file=filename)
}

xgboost_fit<-function(xtrain_dataset){
  library(xgboost)
  set.seed(1)
  gnomad_train<-xgb.DMatrix(data=as.matrix(xtrain_dataset %>% dplyr::select(-outcome)), label=xtrain_dataset$outcome )
  test_ds_tibble<-variants %>% 
    filter(clinvar_holdout_test==TRUE) %>% 
    dplyr::select(any_of(colnames(xtrain_dataset)))
  test_ds<-xgb.DMatrix(data=as.matrix(test_ds_tibble %>% dplyr::select(-outcome)), label=test_ds_tibble$outcome)
  watchlist=list(train=gnomad_train, test=test_ds)
  
  params <- list(max_depth = opt$max_depth_param, 
                 subsample=opt$subsample_param, 
                 eta=opt$eta_param, 
                 lambda=opt$lambda_param, 
                 alpha=0, 
                 gamma=opt$gamma_param, 
                 min_child_weight=opt$min_child_weight_param)
  model1<-xgb.train(
    data=gnomad_train,
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
xgboost_save<-function(modelx, filename){
  xgb.save(model=modelx, fname=filename)
}

if (opt$method_pred=="xgboost"){
  fit_model<-xgboost_fit
  predict_model<-xgboost_predict
  save_model<-xgboost_save
} else if (opt$method_pred=="randomforest") {
  fit_model<-ranger_fit
  predict_model<-ranger_predict
  save_model<-ranger_save
}else {
  fit_model<-extraT_fit
  predict_model<-ranger_predict
  save_model<-ranger_save
}


if (opt$k_fold_cross_val == TRUE){
  set.seed(1)
  
  # splitKFold function from existing_scores_glm_functions.R
  variants<-splitKFold(variants, 5)
    
  glm_cv<-data.frame()
  for (h in unique(variants$kfold_index)) {
    variants<-setTrainTestSet(variants, h)
    
    gnomad_train<-variants %>% filter(gnomad_train)
    
    gnomad_model_Alph<-fit_model(gnomad_train %>% dplyr::select(all_of(colsCorrRemoved)))
    variants$predicted_Alph<-predict_model(variants, gnomad_model_Alph) 
    
    glm_AlphCadd<- glm(outcome ~ . , family=binomial(link='logit'),
                            data=variants %>% filter(clinvar_holdout_test) %>% 
                            dplyr::select(outcome, predicted_Alph, CADD_raw) %>%
                            filter(complete.cases(.)))
    
    variants$glm_AlphCadd <-predict(glm_AlphCadd, variants)#prediction with model = non_h_cv_model_glm, test = : predicted__non_h_cv_glm_h_cv_test
    
    clinvar_holdout_test<-variants %>% filter(clinvar_holdout_test)
    clinvar_interim_test<-variants %>% filter(clinvar_interim_test)
    gnomad_train<-variants %>% filter(gnomad_train)
    
    save_tibble <- rbind(save_tibble, 
                         tibble(auc_Alph_train_gnomAD=roc(gnomad_train$outcome, gnomad_train$predicted_Alph)$auc, 
                                Alph_OOB=ifelse((opt$method_pred=="randomforest"), gnomad_model_Alph$prediction.error, NA),
                                auc_CADD_interim_CV=roc(clinvar_holdout_test$outcome, clinvar_holdout_test$CADD_raw)$auc, 
                                auc_Alph_interim_CV=roc(clinvar_holdout_test$outcome, clinvar_holdout_test$predicted_Alph)$auc,
                                auc_CADD_test_CV=roc(clinvar_interim_test$outcome, clinvar_interim_test$CADD_raw)$auc, 
                                auc_Alph_test_CV=roc(clinvar_interim_test$outcome, clinvar_interim_test$predicted_Alph)$auc, 
                                auc_glm_test_CV=roc(clinvar_interim_test$outcome, clinvar_interim_test$glm_AlphCadd)$auc, 
                                sample_num=h,
                                gnomad_train_nrow=nrow(gnomad_train),
                                interim_cv_nrow=nrow(clinvar_holdout_test),
                                test_cv_nrow=nrow(clinvar_interim_test),
                                condition=opt$prefix ))
    
    
   plot(roc(clinvar_holdout_test$outcome, clinvar_holdout_test$CADD_raw), print.auc = TRUE, col = "red")
   plot(roc(clinvar_holdout_test$outcome, clinvar_holdout_test$predicted_Alph), print.auc = TRUE, 
                               col = "green", print.auc.y = .2, add = TRUE)
    
   plot(roc(clinvar_interim_test$outcome, clinvar_interim_test$CADD_raw), print.auc = TRUE, col = "red", add=FALSE)
   plot(roc(clinvar_interim_test$outcome, clinvar_interim_test$glm_AlphCadd), print.auc = TRUE, 
                           col = "blue", print.auc.y = .2, add = TRUE)
   plot(roc(clinvar_interim_test$outcome, clinvar_interim_test$predicted_Alph), print.auc = TRUE, 
                           col = "green", print.auc.y = .4, add = TRUE)
  }
  
}else if (opt$full_model==TRUE){
  var_full_model<-variants %>% 
    filter(gnomadSet==1)%>%
    dplyr::select(all_of(colsCorrRemoved), "outcome") 

  gnomad_model_Alph<-fit_model(var_full_model)
  var_full_model$predicted_Alph<-predict_model(var_full_model, gnomad_model_Alph)
  variants$predicted_Alph<-predict_model(variants, gnomad_model_Alph)
  
  interim_dataset<-variants %>% 
    filter(clinvar_interim_test)
  
  save_tibble <- tibble(auc_Alph_train_gnomAD=roc(var_full_model$outcome, var_full_model$predicted_Alph)$auc,
                        Alph_OOB=ifelse((opt$method_pred %in% c("randomforest","extratree")), gnomad_model_Alph$prediction.error, NA),
                        gnomad_train_nrow=nrow(var_full_model),
                        gnomad_train_nProteins=train_dataset<-length(unique((variants %>% filter(gnomadSet==1))$Uniprot_acc_split)),
                        auc_CADD_interim_CV=roc(interim_dataset$outcome, interim_dataset$CADD_raw)$auc, 
                        auc_Alph_interim_CV=roc(interim_dataset$outcome, interim_dataset$predicted_Alph)$auc,
                        interim_nrow=nrow(interim_dataset),
                        interim_nProteins=length(unique(interim_dataset$Uniprot_acc_split))
                        )
  
}else{
  train_dataset<-variants %>% 
    filter(gnomad_train)%>%
    dplyr::select(all_of(colsCorrRemoved), "outcome")
  
  gnomad_model_Alph<-fit_model(train_dataset)
  
  train_dataset$predicted_Alph<-predict_model(train_dataset, gnomad_model_Alph)
  variants$predicted_Alph<-predict_model(variants %>% dplyr::select(all_of(colsCorrRemoved), "outcome"), gnomad_model_Alph)
  
  interim_dataset<-variants %>% 
    filter(clinvar_interim_test)
  
  test_dataset<-variants %>% 
    filter(clinvar_holdout_test)
  
  test_dataset_gnomad<-variants %>% 
    filter(gnomad_holdout_test) 
  
  roc_rose <- plot(roc(train_dataset$outcome, train_dataset$predicted_Alph), print.auc = TRUE, col = "blue")
  legend(0.3, 0.3, legend=c("Alph"),
         col=c("blue"), lty=1, cex=0.8)
  title(main = "train set")
  
  roc_rose <- plot(roc(interim_dataset$outcome, interim_dataset$CADD_raw), print.auc = TRUE, col = "red")
  roc_rose <- plot(roc(interim_dataset$outcome, interim_dataset$predicted_Alph), print.auc = TRUE, 
                   col = "blue", print.auc.y = .4, add = TRUE)
  legend(0.3, 0.3, legend=c("CADD", "Alph"),
         col=c("red", "blue"), lty=1, cex=0.8)
  title(main = "interim set")
  
  roc_rose <- plot(roc(test_dataset$outcome, test_dataset$CADD_raw), print.auc = TRUE, col = "red")
  roc_rose <- plot(roc(test_dataset$outcome, test_dataset$predicted_Alph), print.auc = TRUE, 
                   col = "blue", print.auc.y = .4, add = TRUE)
  legend(0.3, 0.3, legend=c("CADD", "Alph"),
         col=c("red", "blue"), lty=1, cex=0.8)
  title(main = "test set")
  
  roc_rose <- plot(roc(test_dataset_gnomad$outcome, test_dataset_gnomad$CADD_raw), print.auc = TRUE, col = "red")
  roc_rose <- plot(roc(test_dataset_gnomad$outcome, test_dataset_gnomad$predicted_Alph), print.auc = TRUE, 
                   col = "blue", print.auc.y = .4, add = TRUE)
  legend(0.3, 0.3, legend=c("CADD", "Alph"),
         col=c("red", "blue"), lty=1, cex=0.8)
  title(main = "test set genes, gnomad variants")
  
  
  save_tibble <- tibble(auc_Alph_train_gnomAD=roc(train_dataset$outcome, train_dataset$predicted_Alph)$auc, 
                        Alph_OOB=ifelse((opt$method_pred!="xgboost"), gnomad_model_Alph$prediction.error, NA),
                        auc_CADD_interim_CV=roc(interim_dataset$outcome, interim_dataset$CADD_raw)$auc, 
                        auc_Alph_interim_CV=roc(interim_dataset$outcome, interim_dataset$predicted_Alph)$auc,
                        auc_CADD_test_CV=roc(test_dataset$outcome, test_dataset$CADD_raw)$auc,
                        auc_Alph_test_CV=roc(test_dataset$outcome, test_dataset$predicted_Alph)$auc, 
                        auc_CADD_test_gnomad=roc(test_dataset_gnomad$outcome, test_dataset_gnomad$CADD_raw)$auc, 
                        auc_Alph_test_gnomad=roc(test_dataset_gnomad$outcome, test_dataset_gnomad$predicted_Alph)$auc, 
                        gnomad_train_nrow=nrow(train_dataset),
                        interim_nrow=nrow(interim_dataset),
                        test_nrow=nrow(test_dataset),
                        test_gnomad_nrow=nrow(test_dataset_gnomad),
                        gnomad_train_nProteins=train_dataset<-length(unique((variants %>% 
                                                                               filter(gnomad_train))$Uniprot_acc_split)),
                        interim_nProteins=length(unique(interim_dataset$Uniprot_acc_split)),
                        test_nProteins=length(unique(test_dataset$Uniprot_acc_split)),
                        test_gnomad_nProteins=length(unique(test_dataset_gnomad$Uniprot_acc_split)),
                        condition=opt$prefix )
}

if (opt$write_dataset){
  write_csv(x=variants, file=paste0(opt$prefix,"_variants.csv.gz"))
}

if (opt$method_pred %in% c("randomforest", "extratree")){
  var_imp<-tibble(importance=as.vector(gnomad_model_Alph$variable.importance), 
                  variable=names(gnomad_model_Alph$variable.importance))
  write_tsv(x=var_imp,
            file = paste0(opt$prefix,"_", opt$importance, "_importance.tsv"))
}

if (opt$write_model){
  save_model(gnomad_model_Alph, file=paste0(opt$prefix,"_written_full_model.RData"))
  saveRDS(colsCorrRemoved, file=paste0(opt$prefix,"_colnames_to_use.RData"))
  saveRDS(toAS_properties,file=paste0(opt$prefix,"_toAS_properties.RData") )
}

write_tsv(x=save_tibble, file=paste0(opt$prefix, "_results.tsv"))


