library(tidyverse)
library(readxl)
library(optparse)

option_list = list(
  make_option(c("-p", "--input_preprocessed"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/SSH/gnomad_extracted_prepro.csv.gz", 
              help="location of csv.gz file with preprocessed, uncalibrated variants"),
  make_option(c("-r", "--input_recalibrated"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/SSH/gnomad_extracted_prepro_recal.csv.gz", 
              help="location of csv.gz file with preprocessed, recalibrated variants"),
  make_option(c("-o", "--out_folder"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/local_files", 
              help="name of folder to store output")

  )

opt = parse_args(OptionParser(option_list=option_list))

variants_pre<-read_csv(opt$input_preprocessed, na=c(".","NA"))
variants_recal<-read_csv(opt$input_recalibrated, na=c(".","NA"))

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

counts_pre<- variants_pre %>%
  filter(!(pure_cv18_to_21_gene==TRUE) & gnomadSet==TRUE)%>%
  group_by(from_AS) %>%
  count(from_AS)

from_AS_pre<-variants_pre %>%
  filter(!(pure_cv18_to_21_gene==TRUE) & gnomadSet==TRUE)%>%
  group_by(from_AS) %>%
  select(from_AS, outcome) %>%
  summarise(frac=mean(outcome))

plot_pre<- ggplot(from_AS_pre,aes(x=reorder(from_AS, -frac), y=frac)) + 
  geom_bar(stat = "identity")+
  geom_text(aes(x=reorder(from_AS, -frac), y=frac, label=counts_pre$n),
            position = position_dodge(width = 1),
            vjust = -0.5, 
            size = 2.5, colour = "black", check_overlap = TRUE)

ggsave("barplot_preprocessed.pdf", width=200, height=100, units="mm", dpi=300)

####recalibrated
counts_recal<- variants_recal %>%
  filter(!(pure_cv18_to_21_gene==TRUE) & gnomadSet==TRUE)%>%
  group_by(from_AS) %>%
  count(from_AS)

from_AS_recal<-variants_recal %>% 
  filter(!(pure_cv18_to_21_gene==TRUE) & gnomadSet==TRUE)%>%
  group_by(from_AS) %>%
  select(from_AS, outcome) %>%
  summarise(frac=mean(outcome))


plot_recal<-ggplot(from_AS_recal,aes(x=factor(from_AS,level=c("TYR","PHE","TRP","LYS","LEU","ASP","HIS","CYS","MET","ILE","GLU","GLN","SER","ASN","GLY","PRO","THR","ALA","VAL","ARG")), y=frac)) + 
  geom_bar(stat = "identity")+
  geom_text(aes(x=reorder(from_AS, -frac), y=frac, label=counts_recal$n),
            position = position_dodge(width = 1),
            vjust = -0.5, 
            size = 2.5, colour = "black", check_overlap = TRUE)

ggsave("barplot_recalibrated.pdf", width=200, height=100, units="mm", dpi=300)