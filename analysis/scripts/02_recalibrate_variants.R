library(tidyverse)

variants_org<-read_csv("data/preprocess/variants_preprocessed.csv.gz")

dir.create("data/recalibrate")
setwd("data/recalibrate")

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

new_variants %>%
  filter(gnomadSet==1)%>%
  group_by(from_AS) %>%
  summarize(mean_outcome=mean(outcome))%>%
  arrange(mean_outcome)

write_csv(new_variants, file="variants_preprocessed_recalibrated.csv.gz")


#### trial strict:

amino_acids<-unique(variants_org$from_AS)
new_variants<-data.frame()
for (amino_acid in amino_acids){
  for (amino_acid_to in unique(variants_org[variants_org$from_AS==amino_acid,]$to_AS)){
  temp_patho<-variants_org %>% 
    filter(outcome==1, from_AS==amino_acid, to_AS==amino_acid_to)%>%
    filter(gnomadSet==TRUE)
  
  temp_benign<-variants_org %>% 
    filter(outcome==0, from_AS==amino_acid, to_AS==amino_acid_to) %>%
    filter(gnomadSet==TRUE)
  
  non_training_variants<-variants_org%>%
    filter(from_AS==amino_acid, to_AS==amino_acid_to) %>%
    filter(gnomadSet==FALSE)
  
  number_ben<-nrow(temp_benign)
  number_patho<-nrow(temp_patho)
  

  if (number_patho>number_ben){
    new_temp_patho<-sample_n(temp_patho,number_ben, replace = FALSE)
    new_temp_ben<-temp_benign
  }else{
    new_temp_patho<-temp_patho
    new_temp_ben<-sample_n(temp_benign,number_patho, replace = FALSE)
  }
  new_variants<-rbind(new_variants, new_temp_ben, new_temp_patho, non_training_variants)
}
}

new_variants %>%
  filter(gnomadSet==1)%>%
  group_by(from_AS) %>%
  summarize(mean_outcome=mean(outcome))%>%
  arrange(mean_outcome)

write_csv(new_variants, file="variants_preprocessed_recalibrated_strict.csv.gz")

