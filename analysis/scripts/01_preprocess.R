library(tidyverse)
library(VennDiagram)
library(data.table)

variants<-fread("../extract_features/data/train_testset1/gnomad_extracted.csv.gz", header=TRUE) 
to_AS_table<-read_tsv("config/to_AS_table.txt")
constraint_scores<-read_tsv("config/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
cv_ids18<-read_tsv("config/clinvar_var_ids_2018.txt", col_names=c("id"), col_types=cols("id"= col_character())) 
cv_ids21<-read_tsv("config/clinvar_var_ids_2021.txt", col_names=c("id"), col_types=cols("id"= col_character()))

dir.create("data/preprocess")
setwd("data/preprocess")


colnames(to_AS_table) <- paste(colnames(to_AS_table), "toAS", sep = "_")

variants<-variants %>% 
  left_join(to_AS_table, by=c("to_AS"="to_AS_toAS"))%>%
  mutate(to_AS=toupper(to_AS))%>%
  mutate(from_AS=RESN)%>%
  mutate(var_id_genomic=paste(`#chr`, `pos(1-based)`, sep=":"))%>%
  mutate(var_id_prot=paste(Uniprot_acc_split, RESI, sep=":"))

head_v<-head(variants)

unip_ids<-tibble(uniprot=unlist(strsplit(variants$Uniprot_acc,";")), 
                 ensg=unlist(strsplit(variants$Ensembl_geneid,";"))) %>% distinct()


constraint_scores<-unip_ids %>% left_join(constraint_scores, by=c("ensg"="gene_id"))%>%
  dplyr::select(uniprot, ensg, gene, oe_lof_upper, pLI)

variants<-variants%>% 
  left_join(constraint_scores, by=c("Uniprot_acc_split"="uniprot"))%>%
  mutate(pLI05=as.integer(pLI>0.5))

cv_temp<-(variants %>% filter(gnomadSet==0, outcome==0))$var_id_genomic
gn_temp<-(variants %>% filter(gnomadSet==1, outcome==0))$var_id_genomic
venn.diagram(x=list(cv_temp, gn_temp),category.names = c("ClinVar" , "gnomAD"), filename = 'ClinVar_vs_gnomAD_benign.png')

cv_temp<-(variants %>% filter(gnomadSet==0))$var_id_genomic
gn_temp<-(variants %>% filter(gnomadSet==1))$var_id_genomic
venn.diagram(x=list(cv_temp, gn_temp),category.names = c("ClinVar" , "gnomAD"), filename = 'ClinVar_vs_gnomAD_all.png')

cv_21prot<-cv_ids21 %>% left_join(variants, by=c("id"="var_id_genomic")) %>% dplyr::select(var_id_prot, id, Uniprot_acc_split) %>% filter(!is.na(Uniprot_acc_split))
cv_18prot<-cv_ids18 %>% left_join(variants, by=c("id"="var_id_genomic")) %>% dplyr::select(var_id_prot, id, Uniprot_acc_split) %>% filter(!is.na(Uniprot_acc_split))
genes_CV_new<-cv_21prot[ ((!cv_21prot$Uniprot_acc_split %in% cv_18prot$Uniprot_acc_split) & (!cv_21prot$id %in% cv_18prot$id)),]
genes_CV_old<-cv_21prot[!cv_21prot$Uniprot_acc_split %in% genes_CV_new$Uniprot_acc_split,]
genes_CV_total<-rbind(genes_CV_new, genes_CV_old)

venn.diagram(x=list(cv_18prot$id, cv_21prot$id, gn_temp),category.names = c("ClinVar2018" , "ClinVar2021", "gnomAD"), filename = 'ClinVar2018vs2021_variants.png')
venn.diagram(x=list(cv_18prot$Uniprot_acc_split, cv_21prot$Uniprot_acc_split),category.names = c("ClinVar2018" , "ClinVar2021"), filename = 'ClinVar2018vs2021_genes.png')

variants$cv21_gene<-(variants$Uniprot_acc_split %in% cv_21prot$Uniprot_acc_split)
variants$pure_cv18_gene<-(variants$Uniprot_acc_split %in% cv_18prot$Uniprot_acc_split)

variant_ids_cv18<-(variants %>% filter(pure_cv18_gene==TRUE))$var_id_genomic
variant_ids_cv21<-(variants %>% filter(cv21_gene==TRUE))$var_id_genomic

variant_ids_gnomAD<-(variants %>% filter(gnomadSet==TRUE))$var_id_genomic

variants$pure_cv18_to_21_gene<-(variants$cv21_gene & (!variants$pure_cv18_gene)) & !(variants$var_id_genomic %in% variant_ids_gnomAD)

variants$gnomad_no_cv21to18<-!variants$pure_cv18_to_21_gene
variants$clinvar_no_cv21to18_no_gnomad<-!variants$pure_cv18_to_21_gene & !(variants$var_id_genomic %in% variant_ids_gnomAD)
variants$cv18_to_21_CV_test<-(variants$pure_cv18_to_21_gene & variants$gnomadSet==0 & !(variants$var_id_genomic %in% variants$var_id_genomic[variants$pure_cv18_to_21_gene==FALSE]))
variants$cv18_to_21_noCV_just_gnomad<-(!variants$cv18_to_21_CV_test & (variants$gnomadSet==1) & !(variants$var_id_genomic %in% variants$var_id_genomic[variants$cv18_to_21_CV_test==TRUE]) & variants$pure_cv18_to_21_gene)

gnomad_no_cv21to18_temp<-(variants %>% filter(gnomad_no_cv21to18==TRUE))$var_id_genomic
clinvar_no_cv21to18_no_gnomad_temp<-(variants %>% filter(clinvar_no_cv21to18_no_gnomad==TRUE))$var_id_genomic
cv18_to_21_CV_test_temp<-(variants %>% filter(cv18_to_21_CV_test==TRUE))$var_id_genomic
cv18_to_21_noCV_just_gnomad_temp<-(variants %>% filter(cv18_to_21_noCV_just_gnomad==TRUE))$var_id_genomic

venn.diagram(x=list(gnomad_no_cv21to18_temp, clinvar_no_cv21to18_no_gnomad_temp, cv18_to_21_CV_test_temp,cv18_to_21_noCV_just_gnomad_temp ),category.names = c("gnomad_no_cv21to18_temp" , "clinvar_no_cv21to18_no_gnomad_temp", "cv18_to_21_CV_test_temp", "cv18_to_21_noCV_just_gnomad_temp"), filename = 'Train_test_sets.png')

colnames(variants)<-gsub("++","..",colnames(variants), fixed=TRUE)
gc()

write_csv(variants, file="variants_preprocessed.csv.gz")

