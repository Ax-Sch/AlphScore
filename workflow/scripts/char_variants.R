library(tidyverse)
library(readxl)
library(optparse)

option_list = list(
  make_option(c("-c", "--csv_location"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/SSH/variants_preprocessed_samp.csv", 
              help="location of csv.gz file"),
  make_option(c("-o", "--out_folder"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/local_files", 
              help="name of folder to store output")
  
)

opt = parse_args(OptionParser(option_list=option_list))

variants<-read_csv(opt$csv_location, na=c(".","NA"))

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

variants_gnomad<- variants %>% filter(gnomadSet==1)
variants_cv<- variants %>% filter(gnomadSet==0)
variants_unique_genes<-unique(variants$var_id_genomic)
variants_unique_proteins<-unique(variants$var_id_prot)
gnomad_outc1<-variants%>%filter(gnomadSet==1 & outcome==1)
gnomad_outc0<-variants%>%filter(gnomadSet==1 & outcome==0)
cv_outc1<-variants%>%filter(gnomadSet==0 & outcome==1)
cv_outc0<-variants%>%filter(gnomadSet==0 & outcome==0)
cv18_to_21_CV<- variants%>%filter(clinvar_holdout_test==1)
gene_no<-unique(variants$gene)
hold_out_genes<- variants%>%filter(hold_out_genes==1)
train1<-variants%>%filter(!(hold_out_genes==TRUE) & gnomadSet==TRUE)
test1<-variants%>%filter(clinvar_interim_test==TRUE)
test2<-variants %>% filter(clinvar_holdout_test==TRUE)
rest_train1<-variants%>%filter((hold_out_genes==TRUE) & gnomadSet==TRUE)
ensemblID_gene<-unique(variants$Ensembl_geneid)
ensemblID_protein<-unique(variants$Ensembl_proteinid)
uniprot<-unique(variants$Uniprot_entry)
uniprot_acc<-unique(variants$Uniprot_acc)
ensemblID_in_genes<-ensemblID_gene %in% variants_unique_genes == TRUE

var_tibble<-tibble(total_var=nrow(variants),
                        total_gnomad=nrow(variants_gnomad),
                        total_cv=nrow(variants_cv),
                        total_genes=length(variants_unique_genes),
                        total_proteins=length(variants_unique_proteins),                        
                        gnomad_outcome1=nrow(gnomad_outc1),
                        gmoad_outcome0=nrow(gnomad_outc0),
                        cv_outcome1=nrow(cv_outc1),
                        cv_outcome0=nrow(cv_outc0),
                        cv18_to_21_CV=nrow(cv18_to_21_CV),
                        gene_no=length(gene_no),
                        train1=nrow(train1),
                        test1=nrow(test1),
                        test2=nrow(test2),
                        rest_train1= nrow(rest_train1),
                        ensemblIDgene=length(ensemblID_gene),
                        uniprot=length(uniprot),
                        uniprot_acc=length(uniprot_acc),
                        ensemblID_in_genes=length(ensemblID_in_genes))
write_tsv(var_tibble,"characteristics_variants.tsv")
