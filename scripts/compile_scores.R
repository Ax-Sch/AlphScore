library(tidyverse)
library(readxl)
list_uniprot_ids<-read_tsv("resources/mapping_genename_uniprot.tsv")

PREDICTED_PROT_SUFFIX="_w_AlphScore_red_FALSE.csv.gz"

sheets<-c("ADRB2","BRCA1","CALM1", "HRAS", "MAPK1", "P53",  "PTEN", "SUMO1", "TPK1", "TPMT", "UBE2I")
scores_conc<-tibble()

for (gene in sheets){
  scores <- read_excel("resources/Source.xlsx", sheet=gene) %>% 
    select(variant, starts_with("DMS"))%>%
    mutate(gene=gene)%>%
    gather(key="DMS", value="DMS_val", -variant, -gene)
  scores_conc=rbind(scores_conc, scores)
}

#MSH2
scores <- read_excel("resources/MSH2_Jia_2020/1-s2.0-S0002929720304390-mmc2.xlsx", sheet="TableS5_AllVariants")%>%
  rename(variant=Variant, DMS=`LOF score`)%>%
  select(variant, starts_with("DMS"))%>%
  mutate(gene="MSH2")%>%
  gather(key="DMS", value="DMS_val", -variant, -gene)
  scores_conc=rbind(scores_conc, scores)


#Abeta
scores <- read_excel("resources/Abeta_Seuma_2020/elife-63364-supp4-v2.xlsx", sheet="1 aa change")%>%
  mutate(variant=paste0(WT_AA, Pos, Mut))%>%
    rename( DMS_nscore=`nscore`)%>%
    select(variant, starts_with("DMS"))%>%
    mutate(gene="Abeta")%>%
    gather(key="DMS", value="DMS_val", -variant, -gene)
  scores_conc=rbind(scores_conc, scores)
  
  
RASH_path="resources/RASH_Bandaru_2017_regulated/elife-27810-supp1-v2.xlsx"
for (gene in excel_sheets(RASH_path)){
  mut_matrix<-read_excel(RASH_path, sheet=gene)
  amino_acid_pos<-as.vector(colnames(mut_matrix))
  amino_acid_ref<-as.vector(unlist(mut_matrix[1,]))
  
  mut_matrix<-as.tibble(t(mut_matrix[2:nrow(mut_matrix),]))
  colnames(mut_matrix)<-as.vector(unlist(mut_matrix[1,]))
  mut_matrix$pos<-amino_acid_pos
  mut_matrix$ref<-amino_acid_ref
  mut_matrix<-mut_matrix[2:nrow(mut_matrix),] %>% 
    gather(key="alt",value="DMS_val", -pos, -ref)%>%
    filter(ref!=alt)%>%
    mutate(variant=paste0(ref, pos, alt))%>%
    mutate(DMS="DMS", 
           gene="RASH")%>%
    select(variant, gene, DMS, DMS_val)
  scores_conc=rbind(scores_conc, mut_matrix)

}


#VKOR1
scores <- read_excel("resources/VKOR1_Chiasson_2020_activity/elife-58026-fig1-data1-v1.xlsx")%>%
  rename( DMS_abundance_score=`abundance_score`)%>%
  rename( DMS_activity_score=`activity_score`)%>%
  select(variant, starts_with("DMS"))%>%
  mutate(gene="VKOR1")%>%
  gather(key="DMS", value="DMS_val", -variant, -gene)
scores_conc=rbind(scores_conc, scores)

scores_conc<-scores_conc %>% 
  filter(!is.na(DMS_val))%>%
  left_join(list_uniprot_ids, by=c("gene"="gene_name"))%>%
  rename(gene_dms=gene)



##### combine vals
library(data.table)

files<-paste0("data/predicted_prots/", list_uniprot_ids$uniprot_id, PREDICTED_PROT_SUFFIX) 

values_joined<-lapply(files, function(x) { #  mc.cores = 1,
  return(fread(x, na.strings = c("NA","."))  ) } )

values_joined<-rbindlist(values_joined)

values_joined<-values_joined%>%
  mutate(variant=paste0(
str_replace_all(RESN_RESI, c( 
 "ALA"="A", "ARG"="R", "ASN"="N", "ASP"="D",
 "CYS"="C", "GLU"="E", "GLN"="Q", "GLY"="G",
 "HIS"="H", "ILE"="I", "LEU"="L", "LYS"="K",
 "MET"="M", "PHE"="F", "PRO"="P", "SER"="S",
 "THR"="T", "TRP"="W", "TYR"="Y", "VAL"="V")),
str_replace_all(to_AS, c( 
  "ALA"="A", "ARG"="R", "ASN"="N", "ASP"="D",
  "CYS"="C", "GLU"="E", "GLN"="Q", "GLY"="G",
  "HIS"="H", "ILE"="I", "LEU"="L", "LYS"="K",
  "MET"="M", "PHE"="F", "PRO"="P", "SER"="S",
  "THR"="T", "TRP"="W", "TYR"="Y", "VAL"="V"))
))

values_joined<-values_joined %>% 
  left_join(scores_conc, by=c("variant"="variant", "Uniprot_acc_split"="uniprot_id"))%>%
  filter(!is.na(DMS_val))%>%
  mutate(DMS_val=as.double(DMS_val))

dir.create("data/validation_set")
setwd("data/validation_set")

fwrite(x=values_joined, file="validation_set_w_AlphScore.csv.gz")


