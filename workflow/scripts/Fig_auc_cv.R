library(tidyverse)
library(readxl)
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("-t", "--tsv_location"), type="character", default="../results/prediction/pre_final_model_k_fold_results.tsv", 
              help="location of .tsv file"),
  make_option(c("-o", "--out_folder"), type="character", default="../results/prediction/plot_k/", 
              help="name of folder to store output")
)

opt = parse_args(OptionParser(option_list=option_list))

output_loop<-read.table(opt$tsv_location, sep="\t",header=TRUE)

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

# rename?
output_loop_num<-select_if(output_loop, is.numeric)


#output_loop_mean<-data.frame()
output_loop_mean<-output_loop %>%
  summarise_if(is.numeric,mean)%>%
  t()

output_loop_sd<-data.frame()
output_loop_sd<-output_loop%>%
  summarise_if(is.numeric,sd)%>%
  t()
  
df<-data.frame(colnames_cv=colnames(output_loop_num),mean=output_loop_mean,sd=output_loop_sd) 

df_sel<-df%>%
  filter(str_detect(colnames_cv, "^auc"))

df_sel$colnames_cv<-sub("auc_", "", df_sel$colnames_cv)

barplot_cv<-ggplot(df_sel, aes(x=reorder(colnames_cv, -mean), y=mean))+
  geom_bar(stat = "identity",width = 0.75, position="dodge")+
  geom_errorbar(aes(ymin=mean - sd,ymax=mean + sd), width=0.4)+
  geom_text(aes(x=reorder(colnames_cv, -mean), y=mean+0.014, label=round(mean, digits=3)), 
            position = position_dodge(0.9),
            vjust = -0.5, 
            size = 5, colour = "black", check_overlap = FALSE)+
  coord_cartesian(ylim=c(0,1))+
  theme_minimal()+
  labs(x="prediction method",y="mean AUC")

barplot_cv

ggsave(filename= "barplot_mean_auc_crossval.pdf",plot=barplot_cv, width=190, height=100, units="mm", dpi=300)
