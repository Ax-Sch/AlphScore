library(tidyverse)
library(pROC)
library(gridExtra)
library(optparse)
set.seed(1)

option_list = list(
  make_option(c("-t", "--test_dataset"), type="character", default="data/prediction/pre_final_model_test_dataset2.csv.gz", 
              help="csv.gz file, test dataset"),
  make_option(c("-v", "--validation_set"), type="character", default="data/validation_set/validation_set_w_AlphScore.csv.gz", 
              help="Validation set"),
  make_option(c("-o", "--out_folder"), type="character", default="data/analyse_score", 
              help="Output folder"))

opt = parse_args(OptionParser(option_list=option_list))

test_dataset2<-read_csv(opt$test_dataset)
validation_dataset<-read_csv(opt$validation_set)

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)
test_dataset2$multi<-!grepl(fixed("F1-"), test_dataset2$pdb_file)

ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted_Alph, fill=(outcome ==0)))+
  facet_wrap(~multi, scales = "free_y")+
  ggtitle("AlphScore values in relation to\nsingle (left) / multi (right) file proteins")+
  xlab("AlphScore")+
  labs(fill = "(l) benign")+
  theme_minimal()

ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted_Alph, fill=(outcome ==0)))+
  facet_wrap(~protein_mean_b_factor>70)+
  ggtitle("AlphScore values in relation to\nmean pLDDT=<70 (left) / >70 (right) file proteins")+
  xlab("AlphScore")+
  labs(fill = "(l) benign")+
  theme_minimal()


ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted_Alph, fill=(outcome ==0)))+
  facet_wrap(~protein_length>600 )+
  ggtitle("AlphScore values in relation to protein length \n =< 600 (left) >600 (right)")+
  xlab("AlphScore")+
  labs(fill = "(l) benign")+
  theme_minimal()

ggplot(test_dataset2)+
  geom_histogram(aes(x=predicted_Alph, fill=(outcome ==0)))+
  facet_wrap(~pLI>0.5 )+
  ggtitle("AlphScore values in relation to protein pLI \n =< 0.5 (left) >0.5 (right)")+
  xlab("AlphScore")+
  labs(fill = "(l) benign")

ggplot(test_dataset2, aes(x=residue_number/protein_length, y=outcome))+
  geom_point(alpha=0.02, size=5)+
  geom_smooth()+
  ggtitle("Pathogenic variants in relation to relative amino acid position")+
  theme_minimal()

ggplot(test_dataset2, aes(x=residue_number/protein_length, y=predicted_Alph))+
  geom_point(alpha=0.5)+
  geom_smooth()+
  ggtitle("AlphScore values in relation to relative position in Protein")+
  theme_minimal()


ggplot(test_dataset2, aes(x=residue_number, y=outcome))+
  geom_jitter(alpha=0.05, size=2, height=0.1)+
  geom_smooth()+
  ggtitle("Pathogenic variants in relation to absolute amino acid position")+
  coord_cartesian(xlim=c(0,3000))+
  theme_minimal()

ggplot(test_dataset2, aes(x=residue_number, y=predicted_Alph))+
  geom_jitter(alpha=0.05, size=2, height=0.1)+
  geom_smooth()+
  ggtitle("AlphaScore in relation to amino acid position")+
  coord_cartesian(xlim=c(0,3000))+
  theme_minimal()

test_dataset2$CADD_plus_ALPH<-(test_dataset2$predicted_Alph+test_dataset2$CADD_raw)
test_dataset2$REVEL_plus_ALPH<-test_dataset2$predicted_Alph+test_dataset2$REVEL_score
test_dataset2$REVEL_plus_CADD<-test_dataset2$CADD_raw+test_dataset2$REVEL_score

rhos<-tibble()
for (aa in unique(test_dataset2$from_AS)){
interim_dataset<-test_dataset2%>% filter(from_AS == aa)
rhos<-rbind(rhos,
  tibble(
  count=nrow(interim_dataset),
  prop_patho=mean(interim_dataset$outcome),
  AA=aa,
  auc_alph=auc(interim_dataset$outcome, interim_dataset$predicted_Alph),
  auc_CADD_plus_ALPH=auc(interim_dataset$outcome, interim_dataset$CADD_plus_ALPH),
  auc_CADD=auc(interim_dataset$outcome,interim_dataset$CADD_raw),
  auc_REVEL=auc(interim_dataset$outcome,interim_dataset$REVEL_score),
  auc_REVEL_plus_ALPH=auc(interim_dataset$outcome,interim_dataset$REVEL_plus_ALPH),
  auc_REVEL_plus_CADD=auc(interim_dataset$outcome,interim_dataset$REVEL_plus_CADD)
  ))
}

#rhos<-rhos %>% 
#  arrange(auc_alph)
#rhos$AA<-factor(rhos$AA, levels = unique(rhos$AA))


allAA_aucs<-tibble(
  scores=c("AlphScore", "CADD", "AlphScore + CADD", "REVEL", "REVEL + AlphScore", "REVEL + CADD"),
  values=c(auc(test_dataset2$outcome, test_dataset2$predicted_Alph),
           auc(test_dataset2$outcome, test_dataset2$CADD_raw),
           auc(test_dataset2$outcome, test_dataset2$CADD_plus_ALPH),
           auc(test_dataset2$outcome, test_dataset2$REVEL_score),
           auc(test_dataset2$outcome, test_dataset2$REVEL_plus_ALPH),
           auc(test_dataset2$outcome, test_dataset2$REVEL_plus_CADD) )
)

p1<-ggplot(rhos)+
  geom_col(aes(x=AA, y=auc_alph))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("AUC AlphScore")+
  xlab("Reference amino acid")+
  geom_hline(yintercept=allAA_aucs$values[1])+
  coord_cartesian(ylim=c(0.5,1))

p2<-ggplot(rhos)+
  geom_col(aes(x=AA, y=auc_CADD))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("AUC CADD")+
  xlab("Reference amino acid")+
  geom_hline(yintercept=allAA_aucs$values[2])+
  coord_cartesian(ylim=c(0.5,1))

p3<-ggplot(rhos)+
  geom_col(aes(x=AA, y=auc_CADD_plus_ALPH))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("AUC AlphScore + CADD")+
  xlab("Reference amino acid")+
  geom_hline(yintercept=allAA_aucs$values[3])+
  coord_cartesian(ylim=c(0.5,1))

p4<-ggplot(rhos)+
  geom_col(aes(x=AA, y=auc_REVEL))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("AUC REVEL")+
  xlab("Reference amino acid")+
  geom_hline(yintercept=allAA_aucs$values[4])+
  coord_cartesian(ylim=c(0.5,1))

p5<-ggplot(rhos)+
  geom_col(aes(x=AA, y=auc_REVEL_plus_ALPH))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("AUC AlphScore + REVEL")+
  xlab("Reference amino acid")+
  geom_hline(yintercept=allAA_aucs$values[5])+
  coord_cartesian(ylim=c(0.5,1))

p6<-ggplot(rhos)+
  geom_col(aes(x=AA, y=auc_REVEL_plus_CADD))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("AUC CADD + REVEL")+
  xlab("Reference amino acid")+
  geom_hline(yintercept=allAA_aucs$values[6])+
  coord_cartesian(ylim=c(0.5,1))

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

ggplot(allAA_aucs)+
  geom_col(aes(x=scores, y=values))+
  coord_cartesian(ylim=c(0.5,1))



validation_dataset<-validation_dataset[is.na(validation_dataset$gnomAD_genomes_AC) & is.na(validation_dataset$gnomAD_genomes_AC),]

ggplot(validation_dataset)+
  geom_point(aes(x=CADD_raw, y=DMS_val), alpha=0.1)+
  facet_wrap(~ paste(DMS, gene_dms), scales="free")

ggplot(validation_dataset)+
  geom_point(aes(x=AlphScore, y=DMS_val), alpha=0.1)+
  facet_wrap(~ paste(DMS, gene_dms), scales="free")


validation_dataset$REVEL_plus_ALPH<-13.5*validation_dataset$AlphScore + 6.3*validation_dataset$REVEL_score
validation_dataset$CADD_plus_ALPH<-18*validation_dataset$AlphScore + 1.3*validation_dataset$CADD_raw
validation_dataset$REVEL_plus_CADD<-5.8*validation_dataset$REVEL_score+0.5*validation_dataset$CADD_raw

validation_dataset$DEOGEN2_score<-str_replace(validation_dataset$DEOGEN2_score, fixed(";."), "")
validation_dataset$DEOGEN2_score<-(str_replace(validation_dataset$DEOGEN2_score, fixed(".;"), ""))

median_score<-function(score_field){
  score_as_list<-str_split(score_field, ";", simplify = TRUE)
  score_as_list[score_as_list=="."]<-NA
  return(median(score_as_list, na.rm=TRUE) )
}

validation_dataset$DEOGEN2_score_med<-sapply(1:nrow(validation_dataset), function(x) { median_score(validation_dataset$DEOGEN2_score[x])})
validation_dataset$DEOGEN2_score_med<-as.numeric(validation_dataset$DEOGEN2_score_med)
validation_dataset$DEOGEN2_plus_ALPH<-validation_dataset$DEOGEN2_score_med+validation_dataset$AlphScore
validation_dataset$DEOGEN2_plus_CADD<-validation_dataset$DEOGEN2_score_med+validation_dataset$CADD_raw
validation_dataset$DEOGEN2_plus_REVEL<-validation_dataset$DEOGEN2_score_med+validation_dataset$REVEL_score
validation_dataset$DEOGEN2_plus_REVEL_plus_ALPH<-validation_dataset$DEOGEN2_score_med+validation_dataset$REVEL_score+validation_dataset$AlphScore


SCORES<-c("AlphScore", "CADD_raw", "CADD_plus_ALPH", "REVEL_score", "REVEL_plus_ALPH", "REVEL_plus_CADD", "DEOGEN2_score_med", "DEOGEN2_plus_ALPH", "DEOGEN2_plus_CADD", "DEOGEN2_plus_REVEL", "DEOGEN2_plus_REVEL_plus_ALPH")
spearmans_joined<-tibble()
for (un_ID in unique(validation_dataset$Uniprot_acc_split)){
  temp_values_joined<-validation_dataset %>% filter(Uniprot_acc_split==un_ID)
  for (dms in unique(temp_values_joined$DMS)){
    temp_inner_loop_values_joined<-temp_values_joined %>% filter(DMS==dms)
    for (score in SCORES){
     score_vector<-as.vector(unlist(temp_inner_loop_values_joined %>% select(one_of(score))))
      spearmans<-tibble(
      gene=unique(temp_inner_loop_values_joined$gene_dms), 
      UP_ID=un_ID, 
      DMS=dms, 
      spearm=ifelse(sum(!is.na(score_vector))>2,
      cor.test(temp_inner_loop_values_joined$DMS_val, score_vector,
               method = "spearman", na.action="na.ommit")$estimate,NA),
      method=score)
        spearmans_joined=rbind(spearmans_joined,spearmans)
      }
  }
}

spearmans_joined<-spearmans_joined %>% mutate(method=factor(method, levels= SCORES))

spearmans_joined_spread<-spearmans_joined%>% spread(method, spearm)

spearmans_joined_summarized<-spearmans_joined %>% 
  group_by(method) %>%
  summarise(mean_abs_correlation=mean(abs(spearm), na.rm = TRUE))%>%
  mutate(method=factor(method, levels= SCORES))

plot_spearmans<-ggplot(spearmans_joined, aes(x=method, y=abs(spearm)))+
  stat_summary(fun.y = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.4)+
  geom_jitter(width=0.15, color="darkblue")+
  xlab("AlphScore")+
  labs(fill = "(l) benign")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1))

plot_spearmans

ggsave(filename= "spearman_plot.pdf", plot=plot_spearmans)

