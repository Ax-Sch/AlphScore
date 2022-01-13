library(tidyverse)
library(pROC)
library(readxl)
library(ranger)
library(optparse)

option_list = list(
  make_option(c("-i", "--input_impurity"), type="character", default="data/prediction_final/final_regular_impurity_importance.tsv", 
              help="location of .tsv file when importance set to impurity"),
  make_option(c("-m", "--input_permutation"), type="character", default="data/prediction_final/final_regular_permutation_importance.tsv", 
              help="location of .tsv file when importance set to impurity"),
  make_option(c("-o", "--out_folder"), type="character", default="data/plot_k/", 
              help="name of folder to store output"),
  make_option(c("-p", "--prefix"), type="character", default="base_tree_2000_10", 
              help="Prefix for output")
  )

opt = parse_args(OptionParser(option_list=option_list))

importance_table_imp<-read.table(opt$input_impurity, sep="\t",header=TRUE)
importance_table_per<-read.table(opt$input_permutation, sep="\t",header=TRUE)

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)


plot_variable_importance<-function(var_imp_ds){
var_importance_plot<-ggplot(var_imp_ds, aes(x=reorder(variable,-importance), y=importance))+ 
  geom_bar(stat="identity", position="dodge")+ 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=10),
        axis.text.y = element_text(color="black", size=10))+
  labs(y = "Variable Importance", size=12)+
  labs(x = "")
return(var_importance_plot)
}

var_importance_imp_plot<-plot_variable_importance(var_imp_imp)

ggsave(filename=paste0(opt$prefix,"_importance_impurity.pdf"), 
       plot=var_importance_imp_plot, width = 6, height = 6)


#Variable importance mode = permutation
var_imp_per<-importance_table_per[order(-importance_table_per$importance),]%>%
  slice(1:25)

var_importance_per_plot<-plot_variable_importance(var_imp_per)

ggsave(filename=paste0(opt$prefix,"_importance_permutation.pdf"), 
       plot=var_importance_per_plot, width = 6, height = 6)
