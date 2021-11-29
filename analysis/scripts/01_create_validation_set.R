library(data.table)
library(tidyverse)
unip_ids<-c("P78352", "P38398", "P51580", "P46937")



files<-paste0("../extract_features/data/combine2_protein/", unip_ids, "_w_AFfeatures.csv.gz") 

values_joined<-lapply(files, function(x) { #  mc.cores = 1,
  return(fread(x, na.strings = c("NA","."))  ) } )

values_joined<-rbindlist(values_joined)

dir.create("data/validation")
setwd("data/validation")



val_files<-c("Brca1.vcf.gz","DLG4.vcf.gz","TPMT.vcf.gz","YAP1.vcf.gz")
base_url<-"https://kircherlab.bihealth.org/download/CADD-development/v1.6/validation/DMS/"


for (val_file in val_files){
download.file(paste0(base_url, val_file), val_file)
}

val_files_joined<-lapply(val_files, function(x) { #  mc.cores = 1,
  return(fread(x, na.strings = c("NA",".") )  ) } )
val_files_joined<-rbindlist(val_files_joined, fill=TRUE)

colnames(val_files_joined)<-c("Chrom","Pos","HGVS_pos","Ref","Alt","Score_exp","Stdev_exp","Stdev2_exp")

validation_set<-values_joined %>% left_join(val_files_joined, by=c("#chr"="Chrom", 
                                                   "hg19_pos(1-based)"="Pos", 
                                                   "ref"="Ref",
                                                   "alt"="Alt"))%>%
  filter(!is.na(Score_exp)) %>%
  mutate(gnomadSet=99,
         outcome=99)


fwrite(validation_set, file="validation_set.csv.gz")

