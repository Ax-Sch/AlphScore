library(tidyverse)
library(gridExtra)
library(optparse)

option_list = list(
  make_option(c("-t", "--variants"), type="character", default="results/prediction_final/pre_final_model_regular_variants.csv.gz", 
              help="csv.gz file, test dataset"),
  make_option(c("-o", "--out_folder"), type="character", default="results/pLDDT", 
              help="Output folder"))

opt = parse_args(OptionParser(option_list=option_list))

variants<-read_csv(opt$variants)

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

gnomAD_vars<-variants %>% 
  filter(gnomadSet==1)%>%
  mutate(outcome_f=ifelse(outcome==1, "proxy pathogenic", "proxy benign"))%>%
  mutate(outcome_f=as.factor(outcome_f))%>%
  mutate(pLDDT_above_70=b_factor>70)

ClinVar_vars<-variants %>% 
  filter(gnomadSet==0)%>%
  mutate(outcome_f=ifelse(outcome==1, "(likely) pathogenic", "(likely) benign"))%>%
  mutate(outcome_f=as.factor(outcome_f))%>%
  mutate(pLDDT_above_70=b_factor>70)


hist_gnomad<-ggplot(gnomAD_vars)+
  geom_histogram(aes(x=b_factor, fill = outcome_f))+
  theme_bw()+
  geom_vline(xintercept=70, color="grey")+
  xlim(15,100)

hist_CV<-ggplot(ClinVar_vars)+
  geom_histogram(aes(x=b_factor, fill = outcome_f))+
  theme_bw()+
  geom_vline(xintercept=70, color="grey")+
  xlim(15,100)

bar_gnomad<-ggplot(gnomAD_vars)+
  geom_bar(aes(fill=outcome_f, x=pLDDT_above_70), position="fill",   width = 0.6)+
  theme_bw()

bar_CV<-ggplot(ClinVar_vars)+
  geom_bar(aes(fill=outcome_f, x=pLDDT_above_70), position="fill", width = 0.6)+
  theme_bw()

arranged_plddt_plots<-grid.arrange(hist_gnomad, hist_CV,bar_gnomad, bar_CV, nrow=2, ncol=2)

ggsave(filename = "arranged_plddt_plots.pdf",  plot=arranged_plddt_plots)

