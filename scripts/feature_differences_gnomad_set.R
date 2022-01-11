library(tidyverse)
library(optparse)
set.seed(1)

option_list = list(
  make_option(c("-v", "--variants"), type="character", default="data/prediction/pre_final_model_variants.csv.gz", 
              help="Validation set"),
  make_option(c("-c", "--colnames_to_use"), type="character", default="data/prediction/pre_final_model_colnames_to_use.RData", 
              help="Output folder"),
  make_option(c("-o", "--out_folder"), type="character", default="data/feature_diff_gnomad", 
              help="Output folder")
)

opt = parse_args(OptionParser(option_list=option_list))

variants<-read_csv(opt$variants)
colnames_to_use<-readRDS(opt$colnames_to_use)

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

gnomad_vars<-variants %>% 
  filter(gnomadSet==1) %>% 
  select(all_of(colnames_to_use), from_AS) %>% 
  select(-(ends_with("toAS") ))

scaled_vars<-as_tibble(scale(gnomad_vars %>% dplyr::select(-outcome, -from_AS) ) )
scaled_vars$outcome<-gnomad_vars$outcome
scaled_vars$from_AS<-gnomad_vars$from_AS

## plot differences
interaction_changes<-scaled_vars %>%
  group_by(from_AS, outcome)%>%
  summarise_all(., mean)

int_changes<-interaction_changes %>% 
  gather(variable, value, -outcome, -from_AS) %>%
  group_by(variable,from_AS) %>%
  summarise(value = (value[outcome == 1] - value[outcome == 0]))

order_vars<-(int_changes %>% 
               dplyr::select(-from_AS) %>% 
               group_by(variable) %>%
               summarise(mean_val=mean(value)) %>% 
               arrange(mean_val))$variable

int_changes <- int_changes %>% 
  mutate(variable=factor(variable, levels=order_vars))

int_changes_matrix<-int_changes %>% spread(variable, value) 
data <- scale(int_changes_matrix %>% dplyr::select(-from_AS))
ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
int_changes_matrix$from_AS[ord]

plot2<-ggplot(int_changes %>% mutate(from_AS=factor(from_AS, levels = int_changes_matrix$from_AS[ord])))+
  geom_tile(aes(x=variable, y=from_AS, fill=value))+
  geom_text(aes(x=variable, y=from_AS, label=round(value,2)), size=1)+
  theme_minimal()+  scale_fill_viridis_c(name="p-pathogenic - \n p-benign")+ #scales::rescale(c(-0.2, 0,0.2, 0.4, 0.6, 0.8,1))
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
  

plot2

ggsave(filename="feature_differences_gnomad_set.pdf", plot=plot2, width=25, height=7)






# plot characteristics of reference amino acids 
interaction_changes<-scaled_vars %>%
  group_by(from_AS)%>%
  summarise_all(., mean)

int_changes<-interaction_changes %>% 
  gather(variable, value, -from_AS) %>%
  group_by(variable,from_AS) 

order_vars<-(int_changes %>% 
               dplyr::select(-from_AS) %>% 
               group_by(variable) %>%
               summarise(mean_val=mean(value)) %>% 
               arrange(mean_val))$variable

int_changes <- int_changes %>% 
  mutate(variable=factor(variable, levels=order_vars))

int_changes_matrix<-int_changes %>% spread(variable, value) 
data <- scale(int_changes_matrix %>% ungroup() %>% dplyr::select(-from_AS))
ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
int_changes_matrix$from_AS[ord]

plot2<-ggplot(int_changes %>% mutate(from_AS=factor(from_AS, levels = int_changes_matrix$from_AS[ord])))+
  geom_tile(aes(x=variable, y=from_AS, fill=value))+
  geom_text(aes(x=variable, y=from_AS, label=round(value,2)), size=1)+
  theme_minimal()+  scale_fill_viridis_c(name="p-pathogenic - \n p-xx")+ #scales::rescale(c(-0.2, 0,0.2, 0.4, 0.6, 0.8,1))
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))

plot2

ggsave(filename="features_of_ref_AAs.pdf", plot=plot2, width=25, height=7)
