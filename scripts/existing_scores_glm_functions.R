# file with functions related to external prediction scores, to be sourced

# function to find the elements in the semicolon separated lists of dbNSFP that correspond to the Alphafold derived score val1;val2;val3
get_index_col<-function(variant_dataset){
  variant_dataset$Uniprot_acc_list<-variant_dataset$Uniprot_acc %>% str_split(";")
  variant_dataset$HGVSp_VEP_list<-variant_dataset$HGVSp_VEP %>% str_split(";")

  variant_dataset$pos_in_Uniprot<-lapply(1:nrow(variant_dataset), function(x) { (variant_dataset$Uniprot_acc_list[[x]] == variant_dataset$Uniprot_acc_split[x])} )
  variant_dataset$pos_in_VEP<-lapply(1:nrow(variant_dataset), function(x) { (variant_dataset$HGVSp_VEP_list[[x]] == variant_dataset$HGVSp_VEP_split[x])}  )
  variant_dataset$pos_in_VEP_and_Uniprot<-lapply(1:nrow(variant_dataset), function(x) {(variant_dataset$pos_in_Uniprot[[x]] & variant_dataset$pos_in_VEP[[x]])}  )
  
  return(variant_dataset$pos_in_VEP_and_Uniprot)
}


# function to retrieve the score that is stored at the positions that corresponds to the Alphafold score
unlist_score<-function(score_col, index_col){
  score_col_split<-str_split(score_col, ";")
  score_splitted<-sapply(1:length(index_col), function(x) {mean(as.numeric(score_col_split[[x]][index_col[[x]]]), na.rm=TRUE )})
  
  return(score_splitted)
}

# Fit multiple GLM models that contain different score-combinations
fit_set_of_models<-function(prefiltered_trainingdataset){
  glm_AlphCadd <- glm(outcome ~ . , family=binomial(link='logit'),
                      data=prefiltered_trainingdataset %>% dplyr::select(outcome, AlphScore, CADD_raw) %>%
                        filter(complete.cases(.)))
  
  glm_AlphRevel <- glm(outcome ~ . , family=binomial(link='logit'),
                       data=prefiltered_trainingdataset %>% dplyr::select(outcome, AlphScore, REVEL_score) %>%
                         filter(complete.cases(.)))
  
  glm_RevelCadd <- glm(outcome ~ . , family=binomial(link='logit'),
                       data=prefiltered_trainingdataset %>% dplyr::select(outcome, CADD_raw, REVEL_score) %>%
                         filter(complete.cases(.)))
  
  glm_AlphRevelCadd <- glm(outcome ~ . , family=binomial(link='logit'),
                           data=prefiltered_trainingdataset %>% dplyr::select(outcome, AlphScore, CADD_raw, REVEL_score) %>%
                             filter(complete.cases(.)))
  
  glm_AlphDeogen <- glm(outcome ~ . , family=binomial(link='logit'),
                        data=prefiltered_trainingdataset %>% dplyr::select(outcome, AlphScore, DEOGEN2_score_med) %>%
                          filter(complete.cases(.)))
  
  glm_CaddDeogen<- glm(outcome ~ . , family=binomial(link='logit'),
                       data=prefiltered_trainingdataset %>% dplyr::select(outcome, CADD_raw, DEOGEN2_score_med) %>%
                         filter(complete.cases(.)))
  
  glm_DeogenRevel <- glm(outcome ~ . , family=binomial(link='logit'),
                         data=prefiltered_trainingdataset %>% dplyr::select(outcome, REVEL_score, DEOGEN2_score_med) %>%
                           filter(complete.cases(.)))
  
  glm_AlphDeogenRevel <- glm(outcome ~ . , family=binomial(link='logit'),
                             data=prefiltered_trainingdataset %>% dplyr::select(outcome, AlphScore, REVEL_score, DEOGEN2_score_med) %>%
                               filter(complete.cases(.)))
  
  glm_AlphCaddDeogen <- glm(outcome ~ . , family=binomial(link='logit'),
                            data=prefiltered_trainingdataset %>% dplyr::select(outcome, AlphScore, CADD_raw, DEOGEN2_score_med) %>%
                              filter(complete.cases(.)))
  
  glm_CaddDeogenRevel <- glm(outcome ~ . , family=binomial(link='logit'),
                             data=prefiltered_trainingdataset %>% dplyr::select(outcome, CADD_raw, DEOGEN2_score_med, REVEL_score) %>%
                               filter(complete.cases(.)))
  
  
  return(list(
    list(glm_AlphCadd, glm_AlphRevel, glm_RevelCadd, glm_AlphRevelCadd, glm_AlphDeogen, 
         glm_CaddDeogen, glm_DeogenRevel, glm_AlphDeogenRevel, glm_AlphCaddDeogen, glm_CaddDeogenRevel),
    
    list("glm_AlphCadd", "glm_AlphRevel", "glm_RevelCadd", "glm_AlphRevelCadd", "glm_AlphDeogen", 
         "glm_CaddDeogen", "glm_DeogenRevel", "glm_AlphDeogenRevel", "glm_AlphCaddDeogen", "glm_CaddDeogenRevel")) 
  )
}

# predict a set of glm models
predict_set_of_models<-function(set_of_models, variants_to_predict){
  # set_of_models[[1]] = models, set_of_models[[2]] = names of models
  for (i in 1:length(set_of_models[[2]])){
    variants_to_predict[,set_of_models[[2]][[i]] ]<-predict(set_of_models[[1]][[i]], variants_to_predict)
  }
  return(variants_to_predict)
}


prepareVariantsForPrediction<-function(varsToUse,to_AS_table_par){
  to_AS_table_mod<-to_AS_table_par
  colnames(to_AS_table_mod) <- paste(colnames(to_AS_table_mod), "toAS", sep = "_")
  
  variants_mod<-varsToUse %>% 
    left_join(to_AS_table_mod, by=c("to_AS"="to_AS_toAS"))%>%
    mutate(to_AS=toupper(to_AS))%>%
    mutate(from_AS=RESN)%>%
    mutate(var_id_genomic=paste(`#chr`, `pos(1-based)`, sep=":"))%>%
    mutate(var_id_prot=paste(Uniprot_acc_split, RESI, sep=":"))
  
  colnames(variants_mod)<-gsub("++","..",colnames(variants_mod), fixed=TRUE)
  return(variants_mod)
}

splitKFold<-function(variants_par, kFold){
  Uniprot_IDs<-tibble(UP_ids=unique(variants_par$Uniprot_acc_split))
  Uniprot_IDs$kfold_index<-as.integer(runif(nrow(Uniprot_IDs),1,kFold+1))
  variants_out<-variants_par %>% 
    left_join(Uniprot_IDs, by=c("Uniprot_acc_split"="UP_ids"), suffix=c("X",""))
  return(variants_out)
}

setTrainTestSet<-function(variants_par, k_fold_num){
  variants_out<-variants_par %>% 
    mutate(train_genes=(kfold_index != k_fold_num))%>%
    mutate(hold_out_genes=(kfold_index == k_fold_num)) %>%
    mutate(gnomad_train=(train_genes & (gnomadSet==1)))
  
  train_var_ids<-variants_out %>% 
    filter(gnomad_train)%>%
    select(var_id_genomic, var_id_prot)
  
  variants_out<-variants_out %>%
    mutate(clinvar_interim_test = train_genes &
          (gnomadSet==0) & 
          !(var_id_prot %in% train_var_ids$var_id_prot) &
          !(var_id_genomic %in% train_var_ids$var_id_genomic) )

  
  interim_var_ids<-variants_out %>% 
    filter(clinvar_interim_test)%>%
    select(var_id_genomic, var_id_prot)
  
  train_and_interim_ids<-rbind(train_var_ids,
                               interim_var_ids)
  
  variants_out<-variants_out %>%
    mutate(clinvar_holdout_test = hold_out_genes &
             (gnomadSet==0) & 
             !(var_id_prot %in% train_and_interim_ids$var_id_prot) &
             !(var_id_genomic %in% train_and_interim_ids$var_id_genomic) )
  
  variants_out<-variants_out %>%
    mutate(gnomad_holdout_test = hold_out_genes &
             (gnomadSet==1) & 
             !(var_id_prot %in% train_and_interim_ids$var_id_prot) &
             !(var_id_genomic %in% train_and_interim_ids$var_id_genomic) )
}
