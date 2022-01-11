library(tidyverse)
library(readxl)
library(gridExtra)
library(optparse)

option_list = list(
  make_option(c("-r", "--input_recalibrated"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/SSH/gnomad_extracted_prepro_rec.csv.gz", 
              help="location of csv.gz file"),
  make_option(c("-e", "--excel_location"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/local_files/available_colnames_regular.xlsx", 
              help="location of excel file"),
  make_option(c("-o", "--out_folder"), type="character", default="C:/Users/Karo.PC-Karo/Master/Semester_3/Lab_rotation_Ludwig/for_Karola/local_files", 
              help="name of folder to store output")
  
)

opt = parse_args(OptionParser(option_list=option_list))

variants<-read_csv(opt$input_recalibrated, na=c(".","NA"))
colnames_usage <- read_excel(opt$excel_location)

dir.create(opt$out_folder, recursive=TRUE)
setwd(opt$out_folder)

write_tsv(as_tibble(colnames(variants)), file="available_colnames_W_surr.tsv")

for_prediction<-(colnames_usage %>% filter(!is.na(for_prediction)))$value

mean_no_na_no_round<-function(x) {(mean(x, na.rm=TRUE))}


mean_properties<-variants %>% filter(b_factor>80) %>% 
  group_by(from_AS, outcome)%>% 
  dplyr::select(for_prediction) %>%
  summarise_all(mean_no_na_no_round)%>%
  ungroup()%>%
  dplyr::select(-starts_with("RESIDUE_NAME"))

mean_properties<-mean_properties %>% 
  gather(variable, value, -outcome, -from_AS) %>%
  group_by(variable,from_AS) %>%
  summarise(value = (value[outcome == 1] - value[outcome == 0]))

####Both Outcomes, total count####
mean_properties_freq<-variants %>%  
  group_by(from_AS, to_AS) %>%
  select(from_AS, to_AS) %>%
  filter(!to_AS=="%3D") %>%
  filter(!to_AS=="???") %>%
  count(from_AS,to_AS,sort=TRUE)

plot_tot_AA<-ggplot(mean_properties_freq, aes(x = from_AS, y = to_AS, fill = n)) +
  geom_tile()+
  theme_minimal()+  
  scale_fill_viridis_c()+
  theme(plot.title= element_text(size = 6),
        axis.text.x = element_text(size=5,angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5, vjust = 0.5, hjust=1),
        legend.key.height = unit(3, 'mm'),
        legend.key.width = unit(2.5, 'mm'),
        legend.title = element_text(size=5), 
        legend.text = element_text(size=5),
        legend.box.spacing = unit(0.1,"mm"),
        axis.title.y=element_text(size=6),
        axis.title.x=element_text(size=6),
        plot.margin=unit(c(1,1,1,1), "mm"))+
  ggtitle("AA exchanges, raw count")+
  geom_text(aes(from_AS, to_AS, label=n), colour = "white", check_overlap = TRUE, size=1.5)

#ggsave("AA_ex_tot.pdf", width=28.575, height=19.05, units="cm", dpi=300)

####Relative fraction outcome1/outcome0####
mean_properties_freq_outcome<-variants %>% 
  filter(!to_AS=="%3D") %>%
  filter(!to_AS=="???") %>%
  group_by(from_AS, to_AS) %>%
  select(from_AS, to_AS, outcome) %>%
  summarise(frac=mean(outcome))

frac_round<- round(mean_properties_freq_outcome$frac, digits = 2)

plot_AA_freq<-ggplot(mean_properties_freq_outcome, aes(x = from_AS, y = to_AS, fill = frac)) +
  geom_tile()+
  theme_minimal()+  
  scale_fill_viridis_c()+
  theme(plot.title= element_text(size = 6),
        axis.text.x = element_text(size=5,angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5, vjust = 0.5, hjust=1),
        legend.key.height = unit(3, 'mm'),
        legend.key.width = unit(2.5, 'mm'),
        legend.title = element_text(size=5), 
        legend.text = element_text(size=5),
        legend.box.spacing = unit(0.1,"mm"),
        axis.title.y=element_text(size=6),
        axis.title.x=element_text(size=6),
        plot.margin=unit(c(1,1,1,1), "mm"))+
  ggtitle("Fraction of AA exchanges with outcome = 1")+
  geom_text(aes(from_AS, to_AS, label=frac_round), colour = "white", check_overlap = TRUE, size=1.5)

#ggsave("AA_ex_frac.pdf",width=28.575, height=19.05, units="cm", dpi=300)

g<-arrangeGrob(plot_tot_AA, plot_AA_freq, ncol=1) #combine both plots in one
ggsave("AA_ex_compared.pdf", width=95, height=120,units="mm", dpi=300,g)
