library(tidyverse)
library(pROC)
library(gridExtra)
library(optparse)
set.seed(1)

option_list = list(
  make_option(c("-v", "--variants"), type="character", default="data/prediction/pre_final_model_variants.csv.gz", 
              help="Validation set"),
  make_option(c("-o", "--out_folder"), type="character", default="data/combine_scores", 
              help="Output folder"))

opt = parse_args(OptionParser(option_list=option_list))

variants<-read_csv(opt$variants)

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)
variants$multi<-!grepl(fixed("F1-"), variants$pdb_file)

#### Analyse score properties in test data set
test_dataset<-variants %>% 
  filter(cv18_to_21_CV_test==TRUE)

p_mult<-ggplot(test_dataset)+
  geom_histogram(aes(x=predicted_Alph, fill=(outcome ==0)))+
  facet_wrap(~multi, scales = "free_y")+
  ggtitle("AlphScore values in relation to\nsingle (left) / multi (right) file proteins")+
  xlab("AlphScore")+
  labs(fill = "(l) benign")+
  theme_minimal()
print(p_mult)
ggsave(filename="p_mult.pdf", p_mult)

p_plDDT<-ggplot(test_dataset)+
  geom_histogram(aes(x=predicted_Alph, fill=(outcome ==0)))+
  facet_wrap(~protein_mean_b_factor>70)+
  ggtitle("AlphScore values in relation to\nmean pLDDT=<70 (left) / >70 (right) file proteins")+
  xlab("AlphScore")+
  labs(fill = "(l) benign")+
  theme_minimal()
print(p_plDDT)
ggsave(filename="p_plDDT.pdf", p_plDDT)

p_length<-ggplot(test_dataset)+
  geom_histogram(aes(x=predicted_Alph, fill=(outcome ==0)))+
  facet_wrap(~protein_length>600 )+
  ggtitle("AlphScore values in relation to protein length \n =< 600 (left) >600 (right)")+
  xlab("AlphScore")+
  labs(fill = "(l) benign")+
  theme_minimal()
print(p_length)
ggsave(filename="p_length.pdf", p_length)

#p_PLI<-ggplot(test_dataset)+
#  geom_histogram(aes(x=predicted_Alph, fill=(outcome ==0)))+
#  facet_wrap(~pLI>0.5 )+
#  ggtitle("AlphScore values in relation to protein pLI \n =< 0.5 (left) >0.5 (right)")+
#  xlab("AlphScore")+
#  labs(fill = "(l) benign")
#print(p_PLI)
#ggsave(filename="p_PLI.pdf", p_PLI)

p_am_pos<-ggplot(test_dataset, aes(x=residue_number/protein_length, y=outcome))+
  geom_point(alpha=0.02, size=5)+
  geom_smooth()+
  ggtitle("Pathogenic variants in relation to relative amino acid position")+
  theme_minimal()
print(p_am_pos)
ggsave(filename="p_am_pos.pdf", p_am_pos)

p_am_pos2<-ggplot(test_dataset, aes(x=residue_number/protein_length, y=predicted_Alph))+
  geom_point(alpha=0.5)+
  geom_smooth()+
  ggtitle("AlphScore values in relation to relative position in Protein")+
  theme_minimal()
print(p_am_pos2)
ggsave(filename="p_am_pos2.pdf", p_am_pos2)

p_resnum_out<-ggplot(test_dataset, aes(x=residue_number, y=outcome))+
  geom_jitter(alpha=0.05, size=2, height=0.1)+
  geom_smooth()+
  ggtitle("Pathogenic variants in relation to absolute amino acid position")+
  coord_cartesian(xlim=c(0,3000))+
  theme_minimal()
print(p_resnum_out)
ggsave(filename="p_resnum_out.pdf", p_resnum_out)

p_resnum_pred<-ggplot(test_dataset, aes(x=residue_number, y=predicted_Alph))+
  geom_jitter(alpha=0.05, size=2, height=0.1)+
  geom_smooth()+
  ggtitle("AlphaScore in relation to amino acid position")+
  coord_cartesian(xlim=c(0,3000))+
  theme_minimal()
print(p_resnum_pred)
ggsave(filename="p_resnum_pred.pdf", p_resnum_pred)



# get DEOGEN2 Scores
variants$DEOGEN2_score<-str_replace(variants$DEOGEN2_score, fixed(";."), "")
variants$DEOGEN2_score<-(str_replace(variants$DEOGEN2_score, fixed(".;"), ""))
median_score<-function(score_field){
  score_as_list<-str_split(score_field, ";", simplify = TRUE)
  score_as_list[score_as_list=="."]<-NA
  score_as_list<-as.double(score_as_list)
  return(median(score_as_list, na.rm=TRUE) )
}
variants$DEOGEN2_score_med<-sapply(1:nrow(variants), function(x) { median_score(variants$DEOGEN2_score[x])})



### combine scores on interim dataset:
train_dataset<-variants %>% 
  filter(train_ds==TRUE)
length(unique(train_dataset$Uniprot_acc_split))

length(unique(variants$Uniprot_acc_split))

interim_test_dataset<-variants %>% 
  filter(CVinterim_no21_18_no_gnomad==TRUE,
         gnomadSet==0)
length(unique(interim_test_dataset$Uniprot_acc_split))



model_glm_AC <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_test_dataset %>% dplyr::select(outcome, predicted_Alph, CADD_raw) %>%
                      filter(complete.cases(.)))

model_glm_AR <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_test_dataset %>% dplyr::select(outcome, predicted_Alph, REVEL_score) %>%
                      filter(complete.cases(.)))

model_glm_RC <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_test_dataset %>% dplyr::select(outcome, CADD_raw, REVEL_score) %>%
                      filter(complete.cases(.)))

model_glm_ARC <- glm(outcome ~ . , family=binomial(link='logit'),
                     data=interim_test_dataset %>% dplyr::select(outcome, predicted_Alph, CADD_raw, REVEL_score) %>%
                       filter(complete.cases(.)))

model_glm_AD <- glm(outcome ~ . , family=binomial(link='logit'),
                    data=interim_test_dataset %>% dplyr::select(outcome, predicted_Alph, DEOGEN2_score_med) %>%
                      filter(complete.cases(.)))


variants$predicted_glm_AC<-predict(model_glm_AC, variants)
variants$predicted_glm_AR<-predict(model_glm_AR, variants)
variants$predicted_glm_RC<-predict(model_glm_RC, variants)
variants$predicted_glm_ARC<-predict(model_glm_ARC, variants)
variants$predicted_glm_AD<-predict(model_glm_AD, variants)

test_dataset<-variants %>% 
  filter(cv18_to_21_CV_test==TRUE)


rhos<-tibble()
for (aa in unique(test_dataset$from_AS)){
  interim_dataset<-test_dataset%>% filter(from_AS == aa)
  rhos<-rbind(rhos,
              tibble(
                prop_patho=mean(interim_dataset$outcome),
                AA=aa,
                auc_alph=auc(interim_dataset$outcome, interim_dataset$predicted_Alph),
                auc_CADD_plus_ALPH=auc(interim_dataset$outcome, interim_dataset$predicted_glm_AC),
                auc_CADD=auc(interim_dataset$outcome,interim_dataset$CADD_raw),
                auc_REVEL=auc(interim_dataset$outcome,interim_dataset$REVEL_score),
                auc_REVEL_plus_ALPH=auc(interim_dataset$outcome,interim_dataset$predicted_glm_AR),
                auc_REVEL_plus_CADD=auc(interim_dataset$outcome,interim_dataset$predicted_glm_RC)
              ))
}

allAA_aucs<-tibble(
  scores=c("AlphScore", "CADD", "AlphScore + CADD", "REVEL", "REVEL + AlphScore", "REVEL + CADD"),
  values=c(auc(test_dataset$outcome, test_dataset$predicted_Alph),
           auc(test_dataset$outcome, test_dataset$CADD_raw),
           auc(test_dataset$outcome, test_dataset$predicted_glm_AC),
           auc(test_dataset$outcome, test_dataset$REVEL_score),
           auc(test_dataset$outcome, test_dataset$predicted_glm_AR),
           auc(test_dataset$outcome, test_dataset$predicted_glm_RC) )
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

combined_plot<-grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)
print(combined_plot)
ggsave(filename="combined_plot.pdf", plot=combined_plot)

aucs<-ggplot(allAA_aucs)+
  geom_col(aes(x=scores, y=values))+
  coord_cartesian(ylim=c(0,1))
print(aucs)
ggsave(filename="aucs.pdf", plot=aucs)


roc_rose <- plot(roc(test_dataset$outcome, test_dataset$CADD_raw), print.auc = TRUE, col = "red")
roc_rose <- plot(roc(test_dataset$outcome, test_dataset$predicted_Alph), print.auc = TRUE, 
                 col = "blue", print.auc.y = .4, add = TRUE)
roc_rose <- plot(roc(test_dataset$outcome, test_dataset$predicted_glm_AC), print.auc = TRUE, 
                 col = "black", print.auc.y = .2, add = TRUE)

roc_rose <- plot(roc(test_dataset$outcome, test_dataset$REVEL_score), print.auc = TRUE, col = "red")
roc_rose <- plot(roc(test_dataset$outcome, test_dataset$predicted_Alph), print.auc = TRUE, 
                 col = "blue", print.auc.y = .4, add = TRUE)
roc_rose <- plot(roc(test_dataset$outcome, test_dataset$predicted_glm_AR), print.auc = TRUE, 
                 col = "black", print.auc.y = .2, add = TRUE)

length(unique(test_dataset$Uniprot_acc_split))

values=c(auc(test_dataset$outcome, test_dataset$predicted_Alph),
         auc(test_dataset$outcome, test_dataset$CADD_raw),
         auc(test_dataset$outcome, test_dataset$predicted_glm_AC),
         auc(test_dataset$outcome, test_dataset$REVEL_score),
         auc(test_dataset$outcome, test_dataset$predicted_glm_AR),
         auc(test_dataset$outcome, test_dataset$predicted_glm_RC) )