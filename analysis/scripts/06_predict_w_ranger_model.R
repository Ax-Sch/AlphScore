library(tidyverse)
library(pROC)
library(readxl)
library(caret)
library(optparse)
library(ranger)
set.seed(1)

option_list = list(
  make_option(c("-c", "--csv_location"), type="character", default="data/validation/validation_set.csv.gz", 
              help="csv.gz file"),
  make_option(c("-m", "--model_location"), type="character", default="data/prediction/base_model_written_full_model.RData", 
              help="Excel file listing columns to use"),
  make_option(c("-o", "--output_file"), type="character", default="data/prediction/predicted.csv.gz", 
              help="Excel file listing columns to use"),
  make_option(c("-u", "--use_cols_file"), type="character", default="data/prediction/colnames_to_use.RData", 
              help="Excel file listing columns to use"),
  make_option(c("-t", "--toAS_properties"), type="character", default="data/prediction/toAS_properties.RData", 
              help="Excel file listing columns to use")
)


opt = parse_args(OptionParser(option_list=option_list))
#DEBUG:
#opt$csv_location="/media/axel/Dateien/Arbeit_Gen/alphafold2/data_from_xcat_v2/variants_preprocessed_recalibrated_v2.csv.gz"

variants<-read_csv(opt$csv_location, na=c(".","NA", NA))
model_to_use<-readRDS(opt$model_location)
colnames_to_use<-readRDS(opt$use_cols_file)
toAS_properties<-readRDS(opt$toAS_properties)

variants<-variants%>%
  left_join(toAS_properties, by=c("from_AS"="from_AS_toAS"))

sel_vars_to<-colnames(toAS_properties)
sel_vars_to<-sel_vars_to[sel_vars_to!="from_AS_toAS"]
variants[, colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])]<-
  variants[, sel_vars_to] - 
  variants[, colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])]

variants$AlphScore<-predict(model_to_use, variants)$predictions

write_csv(x=variants,
          file=opt$output_file)

