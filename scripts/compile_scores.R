library(tidyverse)
library(readxl)
library(data.table)
library(optparse)

option_list = list(
  make_option(c("-m", "--mapping"), type="character", default="resources/mapping_genename_uniprot.tsv", 
              help="location of file that maps gene name in assay to UniProt ID"),
  make_option(c("-p", "--prot_folder"), type="character", default="data/predicted_prots_eval", 
              help="Folder with protein files with AlphScore predictions"),
  make_option(c("-a", "--suffix"), type="character", default="_w_AlphScore_red_FALSE.csv.gz", 
              help="suffix that protein files with AlphScore predictions have")
)
opt = parse_args(OptionParser(option_list=option_list))

list_uniprot_ids<-read_tsv(opt$mapping)
PREDICTED_PROT_SUFFIX=opt$suffix
folder_of_prots=opt$prot_folder

# MAPK1	/ Uniprot P28482 not in mapping file, as this was not mapped between Uniprot and dbNSFP

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
  rename(gene_dms=gene)%>% 
  filter(!is.na(uniprot_id))


files<-paste0(folder_of_prots, list_uniprot_ids$uniprot_id, PREDICTED_PROT_SUFFIX) 

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


