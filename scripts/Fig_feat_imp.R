library(tidyverse)
library(pROC)
library(readxl)
library(ranger)
library(optparse)

option_list = list(
  make_option(c("-i", "--input_impurity"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/local_files/pre_final_model_impurity_importance.tsv", 
              help="location of .tsv file when importance set to impurity"),
  make_option(c("-m", "--input_permutation"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/local_files/pre_final_model_permutation_importance.tsv", 
              help="location of .tsv file when importance set to impurity"),
  make_option(c("-o", "--out_folder"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/local_files", 
              help="name of folder to store output"),
  make_option(c("-p", "--prefix"), type="character", default="base_tree_2000_10", 
              help="Prefix for output")
  )

opt = parse_args(OptionParser(option_list=option_list))

importance_table_imp<-read.table(opt$input_impurity, sep="\t",header=TRUE)
importance_table_per<-read.table(opt$input_permutation, sep="\t",header=TRUE)

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

#Variable importance mode = impurity
var_imp_imp<-importance_table_imp[order(-importance_table_imp$importance),]%>%
  slice(1:25)

var_importance_imp<-ggplot(var_imp_imp, aes(x=reorder(variable,importance), y=importance,fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  ggtitle("Information Value Summary")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")

ggsave(filename=paste0(opt$prefix,"_importance_impurity.pdf"), plot=var_importance_imp, width = 10, height = 10)


#Variable importance mode = permutation
var_imp_per<-importance_table_per[order(-importance_table_per$importance),]%>%
  slice(1:25)

var_importance_per<-ggplot(var_imp_per, aes(x=reorder(variable,importance), y=importance,fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  ggtitle("Information Value Summary")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")

ggsave(filename=paste0(opt$prefix,"_importance_permutation.pdf"), plot=var_importance_per, width = 10, height = 10)

