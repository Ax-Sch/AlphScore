library(tidyverse)
library(pROC)
library(gridExtra)
library(optparse)
library(boot)
source("scripts/existing_scores_glm_functions.R")
BOOT_REPETITIONS=1000
N_CPU=3
COL_ORDER=c("Alph_null","AlphScore","CADD_raw","DEOGEN2_score_med","REVEL_score","glm_AlphCadd","glm_AlphDeogen",
  "glm_AlphRevel","glm_CaddDeogen","glm_DeogenRevel","glm_RevelCadd","glm_AlphDeogenRevel",
  "glm_AlphCaddDeogen","glm_AlphCaddRevel","glm_CaddDeogenRevel")

set.seed(1)


option_list = list(
  make_option(c("-t", "--variants"), type="character", default="results/prediction_final/pre_final_model_regular_variants.csv.gz", 
              help="csv.gz file, test dataset"),
  make_option(c("-v", "--validation_set"), type="character", default="results/validation_set/validation_set_w_AlphScore.csv.gz", 
              help="Validation set"),
  make_option(c("-o", "--out_folder"), type="character", default="results/analyse_score", 
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

variants<-variants%>% 
  mutate(ID=paste(`#chr`, `pos(1-based)`,ref, alt, sep=":"))

### combine scores on interim dataset:
interim_dataset<-variants %>% 
  filter(clinvar_interim_test==TRUE)
gnomad_train_dataset<-variants %>% 
  filter(gnomad_train==TRUE)

validation_dataset<-validation_dataset %>%
  mutate(ID=paste(`#chr`, `pos(1-based)`,ref, alt, sep=":"))%>%
  filter(!ID %in% interim_dataset$ID)%>%
  filter(!ID %in% gnomad_train_dataset$ID)


set_of_models<-fit_set_of_models(interim_dataset)
validation_dataset<-predict_set_of_models(set_of_models, validation_dataset)


selected_scores<-c("HRAS_DMS_g12v",
                   "P53_DMS_WT_Nutlin",
                   "ADRB2_DMS_0.625",          
                   "BRCA1_DMS_E3_high_(a)",    
                   "BRCA1_DMS_(b)",            
                   "MSH2_DMS",      
                   "TPMT_DMS",                 
                   "PTEN_DMS_(a)",             
                   "PTEN_DMS_highqual_(b)",    
                   "SUMO1_DMS_flipped",        
                   "UBE2I_DMS_flipped",       
                   "VKOR1_DMS_abundance_score",
                   "TPK1_DMS_flipped")


validation_dataset<-validation_dataset%>%
  mutate(scoreID=paste(gene_dms, DMS, sep="_"))%>%
  filter(scoreID %in% selected_scores)
  


###### validation Dataset

prot_info_for_join<-validation_dataset %>% 
  distinct(Uniprot_acc_split, protein_length, protein_mean_b_factor)

SCORES<-c("Alph_null","AlphScore", "CADD_raw", "glm_AlphCadd", 
          "REVEL_score", "glm_AlphRevel", "glm_RevelCadd", 
          "DEOGEN2_score_med", "glm_AlphDeogen", "glm_CaddDeogen", 
          "glm_AlphCaddDeogen", "glm_DeogenRevel","glm_AlphDeogenRevel",
          "glm_AlphCaddRevel",    "glm_CaddDeogenRevel")

get_correlations<-function(val_dat){
  

  spearman_joined_private<-tibble()
  
for (un_ID in unique(val_dat$Uniprot_acc_split)){
  temp_values_joined<-val_dat %>% 
    filter(Uniprot_acc_split==un_ID)
  temp_values_joined <- temp_values_joined[complete.cases(temp_values_joined[, SCORES]),]
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
      spearman_joined_private=rbind(spearman_joined_private,spearmans)
      }
  }
}

bad_UP_ids<-(spearman_joined_private %>% 
               filter(number_of_vals<21))$UP_ID

spearman_joined_private<-spearman_joined_private %>% 
  mutate(scoreID=paste(gene, DMS, sep="_"))%>%
  mutate(method=factor(method, levels= SCORES))%>%
  filter(! UP_ID %in% bad_UP_ids)%>%
  left_join(prot_info_for_join, by=c("UP_ID"="Uniprot_acc_split"))%>%
  mutate(abs_spear=abs(spearm)) %>%
  mutate(alphInclude=grepl("Alph", method))

return(spearman_joined_private)
}

spearmans_joined<-get_correlations(validation_dataset)


pairwise.t.test(x = spearmans_joined$abs_spear, g = spearmans_joined$method, paired = TRUE, 
                p.adjust.method = "none", alternative = "greater")
# has been performed first; 
spearmans_joined %>% 
  group_by(gene, UP_ID, DMS) %>% 
  summarise(mean_abs_spear=mean(abs(spearm)))
## selected the scores with the highest overall correlation in mean_corr_per_score

spearmans_joined<-spearmans_joined%>%
  mutate(method=factor(method, levels=COL_ORDER))

plot_spearmans<-ggplot(spearmans_joined, aes(x=method, y=abs(spearm)))+
  stat_summary(fun.y = mean, geom = "bar", width=0.8, color="black", fill="grey") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.4)+
  geom_jitter(width=0.15, color="blue4", alpha=1, size=1)+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1, color="black", size=10),
                          axis.text.y = element_text(color="black", size=10))+
  labs(x = "")+
  labs(y = "absolute Spearman correlation", size=12)

plot_spearmans
ggsave(filename= "spearman_plot.pdf", plot=plot_spearmans, height=5, width=4)

getmeansSD_spearman<-function(spearmans_joined_private){
  table_spearm_private<-spearmans_joined_private %>% 
    mutate(abs_spearman=abs(spearm))%>%
    group_by(method)%>%
    summarise(mean_spearm=mean(abs_spearman), sd=sd(abs_spearman), num=n(), sem=sd(abs_spearman)/sqrt(n())) %>%
    arrange(-mean_spearm)
  return(table_spearm_private)
}

table_spearm<-getmeansSD_spearman(spearmans_joined)
print(table_spearm)

write_tsv(x = table_spearm, file = "table_spearm.tsv")

getMeanAbsCorr<-function(data_priv, i){
  spearman_joined_priv<-get_correlations(data_priv[i,])
  table_spearm_private<-getmeansSD_spearman(spearman_joined_priv)
  return(unlist(table_spearm_private[,"mean_spearm"]))
}

getMeanAbsCorr_scoreids<-function(data_priv, i){
  filtered_data<-tibble()
  
  for (score_id in unlist(data_priv[i,]) ){
  filtered_data<-rbind(filtered_data,
    validation_dataset %>%
      filter(scoreID %in% score_id) )
  }
  
  spearman_joined_priv<-get_correlations(filtered_data)
  table_spearm_private<-getmeansSD_spearman(spearman_joined_priv)
  table_spearm_private$used_scores<-as.numeric(paste(i, collapse = ""))
  # for check print
  #print(data_priv[i,])
  #print(nrow(filtered_data))

  return(table_spearm_private)
}

getMeanAbsCorr_scoreids_return_vector<-function(data_priv, i){
  table_spearm_private<-getMeanAbsCorr_scoreids(data_priv, i)
  correlations_used_scores<-c(unlist(table_spearm_private[,"mean_spearm"]),
                              unlist(table_spearm_private[1,"used_scores"]))

  return(correlations_used_scores)
}

compare_scores<-function(booted_values_private){
  booted_values_private_summarised<-booted_values_private %>%
    rowwise%>%
    select(-starts_with("used_cols_vector"))%>%
    mutate(rowmax=max(across()))%>%
    mutate(glm_AlphDeogenRevel_biggest=glm_AlphDeogenRevel==rowmax)%>%
    mutate(glm_AlphCadd_biggerThanCadd=glm_AlphCadd>CADD_raw)%>%
    mutate(glm_AlphDeogen_biggerThanDeogen=glm_AlphDeogen>DEOGEN2_score_med)%>%
    mutate(glm_AlphRevel_biggerThanRevel=glm_AlphRevel>REVEL_score)
  return(booted_values_private_summarised)
}

get_p_value_table<-function(booted_values_summarised_private){
  relevant_scores<-c("glm_AlphDeogenRevel_biggest","glm_AlphDeogen_biggerThanDeogen","glm_AlphRevel_biggerThanRevel","glm_AlphCadd_biggerThanCadd")
  p_val_tabl<-tibble()
  for (columname in relevant_scores){
    p_val_tabl<- rbind(p_val_tabl, 
                       tibble(name=columname, num_true=sum(unlist(booted_values_summarised_private[,columname])), total_num=nrow(booted_values_summarised_private)))
  }
  return(p_val_tabl)
}



# do bootstrapping for means
mean_original_mat<-getmeansSD_spearman(get_correlations(validation_dataset))
set.seed(1)

boot_meanAbsCor<-boot(data=validation_dataset, statistic=getMeanAbsCorr, R=BOOT_REPETITIONS, parallel="multicore", ncpus=N_CPU)
booted_meanAbsCor_values<-as_tibble(boot_meanAbsCor$t)
colnames(booted_meanAbsCor_values)<-mean_original_mat$method
booted_meanAbsCor_values<-compare_scores(booted_meanAbsCor_values)
pValTable_AbsCor<-get_p_value_table(booted_meanAbsCor_values)
pValTable_AbsCor<-pValTable_AbsCor%>% 
  mutate(pval=1-num_true/total_num)
write_tsv(x=pValTable_AbsCor, "pValTable_AbsCor.tsv")


# bootstrap the choice of scores
score_ids<-as_tibble(unique(validation_dataset$scoreID))
mean_original_mat<-getMeanAbsCorr_scoreids(score_ids,1:13)
set.seed(1)
boot_meanAbsCor<-boot(data=score_ids, statistic=getMeanAbsCorr_scoreids_return_vector, R=BOOT_REPETITIONS, parallel="multicore", ncpus=N_CPU)
booted_meanAbsCor_values<-as_tibble(boot_meanAbsCor$t)
colnames(booted_meanAbsCor_values)<-c(as.character(mean_original_mat$method), "used_cols_vector")
booted_meanAbsCor_values<-compare_scores(booted_meanAbsCor_values)
pValTable_AbsCor<-get_p_value_table(booted_meanAbsCor_values)
pValTable_AbsCor<-pValTable_AbsCor%>% 
  mutate(pval=1-num_true/total_num)
write_tsv(x=pValTable_AbsCor, "pValTable_AbsCor_choice_of_experiments.tsv")
pValTable_AbsCor

rank_score_spearm<-function(spearm_joined_private){
  rank_scored_priv<-spearm_joined_private%>%group_by(scoreID)%>%
    mutate(minCor=min(abs_spear), maxCor=max(abs_spear))%>%
    ungroup()%>%
    mutate(norm_cor=(abs_spear-minCor)/(maxCor-minCor))%>%
    group_by(method,gene)%>%
    summarise(gene_level_cor=mean(norm_cor))%>%
    group_by(method)%>%
    summarise(rank_score=sum(gene_level_cor)/n() , sd_score=sd(gene_level_cor)) %>%
    mutate(sd_score_norm=sd_score/rank_score)%>%
    arrange(-method)
  return(rank_scored_priv)
}

getSpearmanRank<-function(data_priv, i){
  spearman_joined_priv<-get_correlations(data_priv[i,])# %>%
    #filter(method!="Alph_null")
  
  rank_scored<-rank_score_spearm(spearman_joined_priv)
  return(unlist(rank_scored[,2]))
}


# do bootstrapping for rankscores
rank_scored_original_mat<-rank_score_spearm(get_correlations(validation_dataset))
set.seed(1)
boot_rank_rankscore<-boot(data=validation_dataset, statistic=getSpearmanRank, R=BOOT_REPETITIONS, parallel="multicore", ncpus=N_CPU) # parallel="multicore", ncpus=3
booted_rankscores_values<-as_tibble(boot_rank_rankscore$t)
colnames(booted_rankscores_values)<-rank_scored_original_mat$method
booted_rankscores_values<-compare_scores(booted_rankscores_values)
vals<-booted_rankscores_values %>% pivot_longer(cols=1:length(unique(rank_scored_original_mat$method)))

vals<-vals%>%
  mutate(name=factor(name, levels=COL_ORDER))



rankscore_plot<-ggplot(vals)+
  geom_boxplot(aes(x=name,y=value))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=10))
  
ggsave(filename = "rankscore_plot.pdf", plot = rankscore_plot, width=5, height=5.5)

pValTable_RankScore<-get_p_value_table(booted_rankscores_values)
pValTable_RankScore<-pValTable_RankScore%>% 
  mutate(pval=1-num_true/total_num)
pValTable_RankScore
write_tsv(x=pValTable_RankScore, "pValTable_RankScore.tsv")


ggplot(spearmans_joined, aes(x=scoreID, y=abs_spear, color=method, label=method))+
  geom_jitter(position = position_jitter(seed = 1))+
  #geom_text(position = position_jitter(seed = 1))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=10),
        axis.text.y = element_text(color="black", size=10))

