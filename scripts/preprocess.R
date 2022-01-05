library(tidyverse)
library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # prevent VennDiagram to write lots of log messages
library(data.table)
library(optparse)
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
cv_ids18<-read_tsv("resources/clinvar_var_ids_2018.txt", col_names=c("id"), col_types=cols("id"= col_character())) 
cv_ids21<-read_tsv("resources/clinvar_var_ids_2021.txt", col_names=c("id"), col_types=cols("id"= col_character()))

dir.create(dirname(opt$output))
setwd(dirname(opt$output))

colnames(to_AS_table) <- paste(colnames(to_AS_table), "toAS", sep = "_")

variants<-variants %>% 
  left_join(to_AS_table, by=c("to_AS"="to_AS_toAS"))%>%
  mutate(to_AS=toupper(to_AS))%>%
  mutate(from_AS=RESN)%>%
  mutate(var_id_genomic=paste(`#chr`, `pos(1-based)`, sep=":"))%>%
  mutate(var_id_prot=paste(Uniprot_acc_split, RESI, sep=":"))

# check for duplicates
n_occur <- data.frame(table(paste(variants$var_id_genomic, variants$gnomadSet)))
double_ids<-str_split(n_occur[n_occur$Freq > 1,]$Var1, " ", simplify=TRUE)[,1]
double_variants<-variants %>% filter(var_id_genomic %in% double_ids) %>% select(var_id_genomic, Uniprot_acc_split, RESN, RESN_RESI, outcome, gnomadSet, alt)
head(double_variants)
# duplicates are different alternative alleles at one position

cv21s<-cv_ids21 %>% 
  left_join(variants, by=c("id"="var_id_genomic")) %>% 
  dplyr::select(var_id_prot, id, Uniprot_acc_split) %>% 
  filter(!is.na(Uniprot_acc_split))

cv18s<-cv_ids18 %>% 
  left_join(variants, by=c("id"="var_id_genomic")) %>% 
  dplyr::select(var_id_prot, id, Uniprot_acc_split) %>% 
  filter(!is.na(Uniprot_acc_split))

cv_21s_not_18s<-cv21s[ ((!cv21s$Uniprot_acc_split %in% cv18s$Uniprot_acc_split) & (!cv21s$id %in% cv18s$id) & (!cv21s$var_id_prot %in% cv18s$var_id_prot)),]

gnomad_variants<-(variants %>% filter(gnomadSet==0))

venn.diagram(x=list(cv18s$Uniprot_acc_split, cv21s$Uniprot_acc_split, cv_21s_not_18s$Uniprot_acc_split),category.names = c("ClinVar2018\nproteins" , "ClinVar2021\nproteins", "ClinVar2021\nnot 2018\nproteins"), filename = 'ClinVar_2018_2021_proteins.png')
venn.diagram(x=list(cv18s$var_id_prot, cv21s$var_id_prot, cv_21s_not_18s$var_id_prot),category.names = c("ClinVar2018\nvariants" , "ClinVar2021\nvariants", "ClinVar2021\nnot 2018\nvariants"), filename = 'ClinVar_2018_2021_variants.png')
venn.diagram(x=list(cv18s$id, cv21s$id,cv_21s_not_18s$id, gnomad_variants$var_id_genomic),category.names = c("ClinVar2018" , "ClinVar2021","ClinVar2021\nnot 2018", "gnomAD"), filename = 'ClinVar2018vs2021_variants_vs_gnomad.png')

variants$cv21_gene<-(variants$Uniprot_acc_split %in% cv21s$Uniprot_acc_split)
variants$pure_cv18_gene<-(variants$Uniprot_acc_split %in% cv18s$Uniprot_acc_split)
variants$pure_cv18_to_21_gene<-(variants$cv21_gene & (!variants$pure_cv18_gene))




#variants$gnomad_no_cv21to18<-!variants$pure_cv18_to_21_gene
variants$clinvar_no_cv21to18_no_gnomad<-!variants$pure_cv18_to_21_gene & 
  !(variants$var_id_genomic %in% gnomad_variants$var_id_genomic)&
  !(variants$var_id_prot %in% gnomad_variants$var_id_prot)

variants<-variants[!is.na(variants$var_id_genomic),]

variants$cv18_to_21_CV_test<-(variants$pure_cv18_to_21_gene & 
                                variants$gnomadSet==0 & 
                                !(variants$var_id_genomic %in% variants$var_id_genomic[variants$pure_cv18_to_21_gene==FALSE]))


variants$cv18_to_21_noCV_just_gnomad<-((!variants$cv18_to_21_CV_test) & 
                                         variants$pure_cv18_to_21_gene &
                                         (variants$gnomadSet==1) & 
                                         !(variants$var_id_genomic %in% variants$var_id_genomic[variants$cv18_to_21_CV_test==TRUE]) &
                                          !(variants$var_id_prot %in% variants$var_id_prot[variants$cv18_to_21_CV_test==TRUE]) )

variants$train_ds<-!(variants$pure_cv18_to_21_gene==TRUE) & variants$gnomadSet==TRUE & !is.na(variants$gnomadSet)

cv18_to_21_CV_test_temp<-variants[variants$cv18_to_21_CV_test==TRUE,]
clinvar_no_cv21to18_no_gnomad<-variants[variants$clinvar_no_cv21to18_no_gnomad==TRUE,]
cv18_to_21_noCV_just_gnomad_temp<-variants[variants$cv18_to_21_noCV_just_gnomad==TRUE,]
train_ds_temp<-variants[variants$train_ds==TRUE,]
  
venn.diagram(x=list(clinvar_no_cv21to18_no_gnomad$var_id_genomic, cv18_to_21_noCV_just_gnomad_temp$var_id_genomic, cv18_to_21_CV_test_temp$var_id_genomic, train_ds_temp$var_id_genomic ),
             category.names = c("clinvar_no_cv21to18_no_gnomad" , "cv18_to_21_noCV_just_gnomad_temp", "cv18_to_21_CV_test_temp", "potential gnomad train"), filename = 'Train_test_sets.png')
venn.diagram(x=list(clinvar_no_cv21to18_no_gnomad$var_id_prot, cv18_to_21_noCV_just_gnomad_temp$var_id_prot, cv18_to_21_CV_test_temp$var_id_prot, train_ds_temp$var_id_prot),
             category.names = c("clinvar_no_cv21to18_no_gnomad" , "cv18_to_21_noCV_just_gnomad_temp", "cv18_to_21_CV_test_temp", "potential gnomad train"), filename = 'Train_test_sets_prot.png')


colnames(variants)<-gsub("++","..",colnames(variants), fixed=TRUE)
gc()

write_csv(variants, file=basename(opt$output))


