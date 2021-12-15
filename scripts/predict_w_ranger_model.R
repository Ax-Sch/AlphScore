library(readr)
library(dplyr)
library(optparse)
library(ranger)
set.seed(1)

option_list = list(
  make_option(c("-c", "--csv_location"), type="character", default="../extract_features/data/combine2_protein/P38398_w_AFfeatures.csv.gz", 
              help="csv.gz file"),
  make_option(c("-m", "--model_location"), type="character", default="data/prediction/final_written_full_model.RData", 
              help="Excel file listing columns to use"),
  make_option(c("-o", "--output_file"), type="character", default="data/prediction/predicted.csv.gz", 
              help="Excel file listing columns to use"),
  make_option(c("-u", "--use_cols_file"), type="character", default="data/prediction/final_colnames_to_use.RData", 
              help="Excel file listing columns to use"),
  make_option(c("-t", "--toAS_properties"), type="character", default="data/prediction/final_toAS_properties.RData", 
              help="Excel file listing columns to use")
)

opt = parse_args(OptionParser(option_list=option_list))
#DEBUG:
#opt$csv_location="/media/axel/Dateien/Arbeit_Gen/alphafold2/data_from_xcat_v2/variants_preprocessed_recalibrated_v2.csv.gz"

to_AS_table<-read_tsv("resources/to_AS_table.txt")
colnames(to_AS_table) <- paste(colnames(to_AS_table), "toAS", sep = "_")
variants<-read_csv(opt$csv_location, na=c(".","NA", NA))
model_to_use<-readRDS(opt$model_location)
colnames_to_use<-readRDS(opt$use_cols_file)
toAS_properties<-readRDS(opt$toAS_properties)

variants<-variants %>% 
  left_join(to_AS_table, by=c("to_AS"="to_AS_toAS"))%>%
  mutate(to_AS=toupper(to_AS))%>%
  mutate(from_AS=RESN)%>%
  mutate(var_id_genomic=paste(`#chr`, `pos(1-based)`, sep=":"))%>%
  mutate(var_id_prot=paste(Uniprot_acc_split, RESI, sep=":"))

colnames(variants)<-gsub("++","..",colnames(variants), fixed=TRUE)

variants<-variants%>%
  left_join(toAS_properties, by=c("from_AS"="from_AS_toAS"))

colnames_to_use<-colnames_to_use[colnames_to_use!="outcome"]

sel_vars_to<-colnames(toAS_properties)
sel_vars_to<-sel_vars_to[sel_vars_to!="from_AS_toAS"]
variants[, colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])]<-
  variants[, sel_vars_to] - 
  variants[, colnames(toAS_properties[,names(toAS_properties) != "from_AS_toAS"])]

variants<-variants[complete.cases(variants[,colnames_to_use]),]

variants$AlphScore<-predict(model_to_use, variants)$predictions

write_csv(x=variants,
          file=opt$output_file)