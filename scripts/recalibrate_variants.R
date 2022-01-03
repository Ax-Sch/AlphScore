library(tidyverse)
library(optparse)
set.seed(1)

option_list = list(
  make_option(c("-i", "--input"), type="character", default="data/train_testset1/gnomad_extracted_prepro.csv.gz", 
              help="csv.gz file"),
  make_option(c("-o", "--output"), type="character", default="data/preprocess/gnomad_extracted_prepro_rec.csv.gz", 
              help="csv.gz file for output")
)
opt = parse_args(OptionParser(option_list=option_list))

base_file_name=tools::file_path_sans_ext(basename(opt$input))
variants_org<-read_csv(opt$input)

dir.create(dirname(opt$output))
setwd(dirname(opt$output))


variants_org %>% group_by(from_AS) %>%
  summarize(mean_outcome=mean(outcome), anz=sum(outcome))%>%
  arrange(mean_outcome)


amino_acids<-unique(variants_org$from_AS)

new_variants<-data.frame()
for (amino_acid in amino_acids){
  temp_patho<-variants_org %>% 
    filter(outcome==1, from_AS==amino_acid)%>%
    filter(gnomadSet==TRUE)

  temp_benign<-variants_org %>% 
    filter(outcome==0, from_AS==amino_acid) %>%
    filter(gnomadSet==TRUE)

  non_training_variants<-variants_org%>%
    filter(from_AS==amino_acid) %>%
    filter(gnomadSet==FALSE)
  
  number_ben<-nrow(temp_benign)
  number_patho<-nrow(temp_patho)
  sample_number<-ifelse(number_patho>number_ben, number_ben,number_patho)
  
  new_temp_patho<-sample_n(temp_patho,sample_number, replace = FALSE)
  
  new_variants<-rbind(new_variants, temp_benign, new_temp_patho, non_training_variants)
  
}

print(new_variants %>%
  filter(gnomadSet==1)%>%
  group_by(from_AS) %>%
  summarize(mean_outcome=mean(outcome))%>%
  arrange(mean_outcome))

print("total number of variants:")
print(nrow(new_variants))



print("Number of variants in gnomAD set after recalibration:")
print(nrow(new_variants %>%
             filter(gnomadSet==1)))

print("Number of variants in gnomAD set before recalibration:")
print(nrow(variants_org%>%
             filter(gnomadSet==1)))

print("Number of variants in ClinVar set:")
print(nrow(new_variants %>%
             filter(gnomadSet==0)))

write_csv(new_variants, file=basename(opt$output))
