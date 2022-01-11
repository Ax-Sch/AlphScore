library(tidyverse)
library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # prevent VennDiagram to write lots of log messages
library(data.table)
library(optparse)
source("scripts/existing_scores_glm_functions.R")
set.seed(1)

option_list = list(
  make_option(c("-i", "--input"), type="character", default="/media/axel/Dateien/Arbeit_Gen/alphafold2/data_from_xcat_v2/gnomad_extracted_v2.csv.gz", 
              help="csv.gz file"),
  make_option(c("-o", "--output"), type="character", default="data/train_testset1/gnomad_extracted_prepro.csv.gz", 
              help="csv.gz file for output")
)
opt = parse_args(OptionParser(option_list=option_list))

variants<-read_csv(opt$input)
to_AS_table<-read_tsv("resources/to_AS_table.txt")

dir.create(dirname(opt$output))
setwd(dirname(opt$output))

# preprocess variants (see sourced file)
variants<-prepareVariantsForPrediction(variants, to_AS_table)
variants<-splitKFold(variants, 5)
variants<-setTrainTestSet(variants, 1)

# check for overlaps
gnomad_train<-variants %>% filter(gnomad_train)
clinvar_holdout_test<-variants %>% filter(clinvar_holdout_test)
clinvar_interim_test<-variants %>% filter(clinvar_interim_test)

venn.diagram(x=list(gnomad_train$Uniprot_acc_split, clinvar_holdout_test$Uniprot_acc_split, clinvar_interim_test$Uniprot_acc_split),category.names = c("gnomAD train" , "ClinVar interim", "ClinVar hold out"), filename = 'Protein_id_level.png')
venn.diagram(x=list(gnomad_train$var_id_prot, clinvar_holdout_test$var_id_prot, clinvar_interim_test$var_id_prot),category.names = c("gnomAD train" , "ClinVar interim", "ClinVar hold out"), filename = 'Protein_variant_id_level.png')
venn.diagram(x=list(gnomad_train$var_id_genomic, clinvar_holdout_test$var_id_genomic, clinvar_interim_test$var_id_genomic),category.names = c("gnomAD train" , "ClinVar interim", "ClinVar hold out"), filename = 'Variant_id_level.png')

# check for duplicates
n_occur <- data.frame(table(paste(variants$var_id_genomic, variants$gnomadSet)))
double_ids<-str_split(n_occur[n_occur$Freq > 1,]$Var1, " ", simplify=TRUE)[,1]
double_variants<-variants %>% filter(var_id_genomic %in% double_ids) %>% select(var_id_genomic, Uniprot_acc_split, RESN, RESN_RESI, outcome, gnomadSet, alt)
head(double_variants)
# duplicates are different alternative alleles at one position

write_csv(variants, file=basename(opt$output))
