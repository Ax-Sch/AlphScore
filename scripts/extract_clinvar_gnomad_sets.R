library(data.table)

files <- list.files(path=snakemake@params[["loaction_csvgzs"]], pattern=".csv.gz", full.names=TRUE, recursive=FALSE)
for (x_files in 1:length(files)){
  if (file.info(files[x_files])$size < 1000){
    files[x_files]<-NA
  }
}
print(files)
files<-files[!is.na(files)]

#testing
#files<-head(files, n=5)

values_joined<-lapply(files, function(x) { #  mc.cores = 1,
  training_data<-fread(x, na.strings = c("NA","."))[!is.na(CADD_raw_sur), ] # load data
  
  singletons <- training_data[( ((is.na(gnomAD_genomes_AC) | gnomAD_genomes_AC<2) & 
                                  gnomAD_exomes_AC==1 & 
                                  gnomAD_exomes_NFE_AC==1) & 
                                 (is.na(`1000Gp3_AC`) | `1000Gp3_AC`==0)  & 
                                 (is.na(ESP6500_AA_AC) | ESP6500_AA_AC==0) & 
                                 (is.na(ESP6500_EA_AC) | ESP6500_AA_AC==0) 
                                )]
  
  frequents<- training_data[gnomAD_exomes_AF>0.001 | 
                              gnomAD_genomes_AF> 0.001]
  
  clinvars<- training_data[clinvar_clnsig %in% c("Likely_pathogenic","Pathogenic","Benign", "Likely_benign")]
  
  return(rbindlist(list(singletons[,`:=`(outcome=1, gnomadSet=1)],
                        frequents[,`:=` (outcome=0, gnomadSet=1)], 
                        clinvars[,`:=` (outcome=ifelse(clinvar_clnsig %in% c("Likely_pathogenic","Pathogenic"),1,0), gnomadSet=0)]), 
                   fill=TRUE, use.names=TRUE))
})

values_joined<-rbindlist(values_joined)

fwrite(values_joined, file=snakemake@output[["csv_gz_out"]])