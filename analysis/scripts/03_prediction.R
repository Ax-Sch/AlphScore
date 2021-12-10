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
  make_option(c("-p", "--prefix"), type="character", default="base_model", 
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
              help="randomforest: min.node.size, xgboost: min node parameter, extratree: nodesize"),
  make_option(c("-r", "--cor_param"), type="double", default=0.999999, 
              help="Just keep columns that have a lower correlation than this value, general"),
  make_option(c("-b", "--b_factor_param"), type="integer", default=90, 
              help="pLDDT value to filter positions that will be used to generate the average values of amino acids, general"),
  make_option(c("-d", "--min_child_weight_param"), type="integer", default=0, 
              help="min child weight, xgboost"),
  make_option(c("-g", "--gamma_param"), type="integer", default=0, 
              help="gamma value, xgboost"),
  make_option(c("-o", "--out_folder"), type="character", default="data/prediction", 
              help="name of folder to store output"),
  make_option(c("-w", "--write_dataset"), type="logical", default=TRUE, 
              help="Write predictions of test-dataset to file"),
  make_option(c("-v", "--validation_set"), type="character", default="data/preprocess/validation_set.csv.gzpreprocessed.csv.gz", 
              help="validation set to calculate scores of"),
  make_option(c("-k", "--k_fold_cross_val"), type="logical", default=FALSE,
              help="activate k-fold cross validation of the training and test data set")
  
)

opt = parse_args(OptionParser(option_list=option_list))
#DEBUG:
#opt$csv_location="/media/axel/Dateien/Arbeit_Gen/alphafold2/data_from_xcat_v2/variants_preprocessed_recalibrated_v2.csv.gz"


save_tibble<-data.frame()

variants<-read_csv(opt$csv_location, na=c(".","NA", NA))
colnames_usage <- read_excel(opt$excel_locataion, col_types="text")

if (opt$validation_set != ""){
  validation_set<-read_csv(opt$validation_set)
}


dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)
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


ranger_fit<-function(xtrain_dataset){
  library(ranger)
  model1<-ranger(outcome ~ ., 
                 data=xtrain_dataset, 
                 importance="impurity", 
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
                 importance="impurity", 
                 splitrule="extratrees",
                 max.depth=opt$max_depth_param,
                 num.trees = opt$num_trees_param,
                 min.node.size = opt$min_node_param
  ) 
  return(model1)
}

xgboost_fit<-function(xtrain_dataset){
  library(xgboost)
  set.seed(1)
  n_ds<-nrow(xtrain_dataset)
  train_i<-sample(1:n_ds,as.integer(n_ds*0.8))
  train_ds<-xgb.DMatrix(data=as.matrix(xtrain_dataset[train_i,] %>% dplyr::select(-outcome)), label=xtrain_dataset$outcome[train_i] )
  test_ds<-xgb.DMatrix(data=as.matrix(xtrain_dataset[-train_i,] %>% dplyr::select(-outcome)), label=xtrain_dataset$outcome[-train_i] )
  watchlist=list(train=train_ds, test=test_ds)
  
  params <- list(max_depth = opt$max_depth_param, 
                 subsample=opt$subsample_param, 
                 eta=opt$eta_param, 
                 lambda=opt$lambda_param, 
                 alpha=0, 
                 gamma=opt$gamma_param, 
                 min_child_weight=opt$min_child_weight_param)
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
  predict_model<-ranger_predict
}


if (opt$k_fold_cross_val == TRUE){
  set.seed(1)
  genes<-unique(variants$gene)
  genes<-as_tibble(genes)
  genes$index<-as.integer(runif(nrow(genes),1,6))
  
  
  glm_cv<-data.frame()
  
  cc<-complete.cases(variants[, c(colnames_new, "outcome")])
  variants_filtered<-variants[cc, ]
  
  for (h in unique(genes$index)) {
    non_h_gnomad_train<-variants_filtered %>% filter(!gene %in% genes[genes$index==h,]$value) %>% #non_h = every gene with index = non-h (if h=1 then non_h=2,3,4,5 etc.); gnomadSet= 1 --> gnomad as training set
      filter(gnomadSet==1) 
    non_h_cv_test <- variants_filtered %>% filter(!gene %in% genes[genes$index==h,]$value) %>% #non_h = every gene with index = non-h (if h=1 then 2,3,4,5 etc.); gnomadSet= 0 --> cv (ClinVar) as test set
      filter(gnomadSet == 0,
             !var_id_genomic %in% non_h_gnomad_train$var_id_genomic)%>%
      filter(!var_id_prot %in% non_h_gnomad_train$var_id_prot,
             !var_id_genomic %in% non_h_gnomad_train$var_id_genomic)
    
    h_cv_test<-variants_filtered %>% filter(gene %in% genes[genes$index==h,]$value) %>% #h = every gene with index = current h (if h=1 then 1); gnomadSet= 0 --> cv as test set
      filter(gnomadSet == 0)%>%
      filter(!var_id_prot %in% non_h_gnomad_train$var_id_prot,
             !var_id_genomic %in% non_h_gnomad_train$var_id_genomic,
             !var_id_prot %in% non_h_cv_test$var_id_prot,
             !var_id_genomic %in% non_h_cv_test$var_id_genomic)
    
    
    gnomad_model_Alph<-fit_model(non_h_gnomad_train %>% dplyr::select(all_of(colnames_new)))
    non_h_cv_test$predicted_Alph<-predict_model(non_h_cv_test %>% dplyr::select(all_of(colnames_new)), gnomad_model_Alph) #prediction with model = non_h_gnomad_Alph and test = non_h_cv_test: predicted_non_h_gnomad_Alph_non_h_cv_test
    non_h_gnomad_train$predicted_Alph<-predict_model(non_h_gnomad_train %>% dplyr::select(all_of(colnames_new)), gnomad_model_Alph) 
    
    non_h_cv_model_glm<- glm(outcome ~ . , family=binomial(link='logit'), #glm model with non_h_cv, inclusion of CADD score; predicted_non_h_gnomad_Alph_non_h_cv_test
                             data=non_h_cv_test %>% dplyr::select(outcome, predicted_Alph, CADD_raw) %>%
                               filter(complete.cases(.)))
    
    h_cv_test$predicted_Alph <- predict_model(h_cv_test %>% dplyr::select(all_of(colnames_new)), gnomad_model_Alph) #prediction with model = non_h_gnomad_Alph and test = h_cv_test: predicted_non_h_gnomad_Alph_h_cv_test
    h_cv_test$predicted_glm <-predict(non_h_cv_model_glm, h_cv_test)#prediction with model = non_h_cv_model_glm, test = : predicted__non_h_cv_glm_h_cv_test
    
    
    save_tibble <- rbind(save_tibble, 
                         tibble(auc_Alph_non_h_gnomAD=roc(non_h_gnomad_train$outcome, non_h_gnomad_train$predicted_Alph)$auc, 
                                Alph_OOB=ifelse((opt$method_pred=="randomforest"), gnomad_model_Alph$prediction.error, NA),
                                auc_CADD_non_h_CV=roc(non_h_cv_test$outcome, non_h_cv_test$CADD_raw)$auc, 
                                auc_Alph_non_h_CV=roc(non_h_cv_test$outcome, non_h_cv_test$predicted_Alph)$auc,
                                auc_CADD_h_CV=roc(h_cv_test$outcome, h_cv_test$CADD_raw)$auc, 
                                auc_Alph_h_CV=roc(h_cv_test$outcome, h_cv_test$predicted_Alph)$auc, 
                                auc_glm_h_CV=roc(h_cv_test$outcome, h_cv_test$predicted_glm)$auc, 
                                sample_num=h,
                                non_h_gnomad_train_nrow=nrow(non_h_gnomad_train),
                                non_h_cv_test_nrow=nrow(non_h_cv_test),
                                h_cv_test_nrow=nrow(h_cv_test),
                                condition=opt$prefix, 
                                params=I(list(opt))))
    
    
    plot_auc_CADD_non_h <- plot(roc(non_h_cv_test$outcome, non_h_cv_test$CADD_raw), print.auc = TRUE, col = "red")
    plot_auc_Alph_non_h<- plot(roc(non_h_cv_test$outcome, non_h_cv_test$predicted_Alph), print.auc = TRUE, 
                               col = "green", print.auc.y = .2, add = TRUE)
    
    plot_auc_CADD_h <- plot(roc(h_cv_test$outcome, h_cv_test$CADD_raw), print.auc = TRUE, col = "red", add=FALSE)
    plot_auc_Alph_h<- plot(roc(h_cv_test$outcome, h_cv_test$predicted_glm), print.auc = TRUE, 
                           col = "blue", print.auc.y = .2, add = TRUE)
    plot_auc_glm_h <- plot(roc(h_cv_test$outcome, h_cv_test$predicted_Alph), print.auc = TRUE, 
                           col = "green", print.auc.y = .4, add = TRUE)
  }
  
}else{
  train_dataset<-train_dataset[,colnames_new] 
  train_dataset<-train_dataset[complete.cases(train_dataset),]
  
  interim_dataset<-variants %>% 
    filter(clinvar_no_cv21to18_no_gnomad==TRUE)
  interim_dataset<-interim_dataset[complete.cases(interim_dataset[,c(colnames_new, "outcome", "CADD_raw")]), ]
  
  gnomad_model_Alph<-fit_model(train_dataset %>% dplyr::select(all_of(colnames_new)))
  
  ## OOB error:
  #gnomad_model_Alph$prediction.error
  
  train_dataset$predicted_Alph<-predict_model(train_dataset %>% dplyr::select(all_of(colnames_new)), gnomad_model_Alph)
  roc_rose <- plot(roc(train_dataset$outcome, train_dataset$predicted_Alph), print.auc = TRUE, col = "red")
  
  interim_dataset$predicted_Alph<-predict_model(interim_dataset %>% dplyr::select(all_of(colnames_new)), gnomad_model_Alph)
  
  roc_rose <- plot(roc(interim_dataset$outcome, interim_dataset$CADD_raw), print.auc = TRUE, col = "red")
  roc_rose <- plot(roc(interim_dataset$outcome, interim_dataset$predicted_Alph), print.auc = TRUE, 
                   col = "green", print.auc.y = .4, add = TRUE)
  
  
  model_glm <- glm(outcome ~ . , family=binomial(link='logit'),
                data=interim_dataset %>% dplyr::select(outcome, predicted_Alph, CADD_raw) %>%
                  filter(complete.cases(.)))
  
  test_dataset<-variants %>% 
    filter(cv18_to_21_CV_test==TRUE)
  
  test_dataset<-test_dataset[complete.cases(test_dataset[,c(colnames_new, "outcome")]), ]
  test_dataset$predicted_Alph<-predict_model(test_dataset %>% dplyr::select(all_of(colnames_new)), gnomad_model_Alph)
  test_dataset$predicted_scaled<-scale(test_dataset$predicted_Alph)
  
  test_dataset$predicted_glm<-predict(model_glm, test_dataset)
  
  
  roc_rose <- plot(roc(test_dataset$outcome, test_dataset$CADD_raw), print.auc = TRUE, col = "red")
  roc_rose <- plot(roc(test_dataset$outcome, test_dataset$predicted_glm), print.auc = TRUE, 
                   col = "blue", print.auc.y = .2, add = TRUE)
  roc_rose <- plot(roc(test_dataset$outcome, test_dataset$predicted_Alph), print.auc = TRUE, 
                   col = "green", print.auc.y = .4, add = TRUE)
  
  save_tibble <- rbind(save_tibble, 
                       tibble(auc_Alph_train_gnomAD=roc(train_dataset$outcome, train_dataset$predicted_Alph)$auc, 
                              Alph_OOB=ifelse((opt$method_pred=="randomforest"), gnomad_model_Alph$prediction.error, NA),
                              auc_CADD_interim_CV=roc(interim_dataset$outcome, interim_dataset$CADD_raw)$auc, 
                              auc_Alph_interim_CV=roc(interim_dataset$outcome, interim_dataset$predicted_Alph)$auc,
                              auc_CADD_test_CV=roc(test_dataset$outcome, test_dataset$CADD_raw)$auc, 
                              auc_Alph_test_CV=roc(test_dataset$outcome, test_dataset$predicted_Alph)$auc, 
                              auc_glm_test_CV=roc(test_dataset$outcome, test_dataset$predicted_glm)$auc, 
                              gnomad_train_nrow=nrow(train_dataset),
                              interim_nrow=nrow(interim_dataset),
                              test_nrow=nrow(test_dataset),
                              condition=opt$prefix, 
                              params=I(list(opt))))
  
  write.csv2(x=save_tibble, file=paste0(opt$prefix, "_results.tsv"))
  
  ggplot(test_dataset)+
    geom_point(aes(x=predicted_Alph, y=CADD_raw), alpha=0.3)+
    facet_wrap(~outcome)
  
  if (opt$write_dataset){
    write_csv(x=test_dataset, file=paste0(opt$prefix,"_test_dataset2.csv.gz"))
  }
}


if (opt$method_pred=="randomforest"){
  ##### CHECK properties of models
  var_imp<-tibble(importance=as.vector(gnomad_model_Alph$variable.importance), variable=names(gnomad_model_Alph$variable.importance))
  
  var_importance<-ggplot(var_imp, aes(x=reorder(variable,importance), y=importance,fill=importance))+ 
    geom_bar(stat="identity", position="dodge")+ coord_flip()+
    ylab("Variable Importance")+
    xlab("")+
    ggtitle("Information Value Summary")+
    guides(fill=F)+
    scale_fill_gradient(low="red", high="blue")
  
  ggsave(filename=paste0(opt$prefix,"_importance.pdf"), plot=var_importance, height=49)
}



if (opt$validation_set != ""){
  validation_set<-validation_set%>%
    left_join(toAS_properties, by=c("from_AS"="from_AS_toAS"))
  
  validation_set[, colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])]<-
    validation_set[, sel_vars_to] - 
    validation_set[, colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])]
  
  validation_set<-validation_set[complete.cases(validation_set[,c(colnames_new)]), ]
  
  validation_set$predicted_Alph<-predict_model(validation_set %>% dplyr::select(all_of(colnames_new[colnames_new!="outcome"])), gnomad_model_Alph)
  validation_set$predicted_glm<-predict(model_glm, validation_set)
  
  write_csv(x=validation_set, file=paste0(opt$prefix,"_validation_set.csv.gz"))
}

