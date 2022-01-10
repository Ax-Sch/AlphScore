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

from_AS_pre<-variants_pre %>%
  filter(!(pure_cv18_to_21_gene==TRUE) & gnomadSet==TRUE)%>%
  group_by(from_AS) %>%
  summarise(frac=mean(outcome), num=n())

plot_pre<- ggplot(from_AS_pre,aes(x=reorder(from_AS, -frac), y=frac)) + 
  geom_bar(stat = "identity")+
  geom_text(aes(x=reorder(from_AS, -frac), y=frac, label=num),
            position = position_dodge(width = 1),
            vjust = -0.5, 
            size = 2.5, colour = "black", check_overlap = TRUE)+
  theme_minimal()+
  labs(x="reference AA",y="proportion of proxy-pathogenic AAs")

ggsave("barplot_preprocessed.pdf", width=200, height=100, units="mm", dpi=300)

order_AA_pre<- from_AS_pre[order(-from_AS_pre$frac),]
order_AA_pre<-order_AA_pre$from_AS

####recalibrated

from_AS_recal<-variants_recal %>% 
  filter(!(pure_cv18_to_21_gene==TRUE) & gnomadSet==TRUE)%>%
  group_by(from_AS) %>%
  summarise(frac=mean(outcome), num=n())

from_AS_recal<-from_AS_recal[match(order_AA_pre, from_AS_recal$from_AS),]
from_AS_recal$from_AS<-factor(from_AS_recal$from_AS, levels=from_AS_recal$from_AS)

plot_recal<-ggplot(from_AS_recal,aes(x=from_AS, y=frac)) + 
  geom_bar(stat = "identity")+
  geom_text(aes(x=from_AS, y=frac, label=num),
            position = position_dodge(width = 1),
            vjust = -0.5, 
            size = 2.5, colour = "black", check_overlap = TRUE)+
  theme_minimal()+
  labs(x="reference AA",y="proportion of proxy-pathogenic AAs")

# get stats of training / testing data set:


ggsave("barplot_recalibrated.pdf", width=200, height=100, units="mm", dpi=300)
