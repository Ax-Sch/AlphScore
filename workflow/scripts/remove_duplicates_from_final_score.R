library(tidyverse)
library(optparse)

option_list = list(
  make_option(c("-v", "--dup_vars"), type="character", default="results/merge_final/all_possible_values_dups.csv.gz", 
              help="list of duplicated variants, final set"),
  make_option(c("-o", "--out_file"), type="character", default="results/merge_final/deduplicated_vars_bind.tsv", 
              help="Output file name"))

opt = parse_args(OptionParser(option_list=option_list))


duplicted_vars<-read_tsv(opt$dup_vars)

collapse_as_text<-c("Uniprot_acc_split", "HGVSp_VEP_split", "b_factor", "SOLVENT_ACCESSIBILITY_core","in_gnomad_train","in_clinvar_ds")
scoreVals_to_combine<-c("CADD_raw", "REVEL_score", "AlphScore", "glm_AlphCadd", "glm_AlphRevel", "glm_RevelCadd", "glm_AlphRevelCadd","glm_AlphDeogen", "glm_CaddDeogen", "glm_DeogenRevel", "glm_AlphDeogenRevel", "glm_AlphCaddDeogen", "glm_CaddDeogenRevel")
any_action_required<-c(collapse_as_text, scoreVals_to_combine)

cols_vars<-colnames(duplicted_vars)
constant_cols<-cols_vars[!cols_vars %in% any_action_required]

deduplicate_vars<-function(varID_of_duplicate){
  duplicated_rows<-duplicted_vars %>% 
    mutate(temp_var_ID=paste(ID, genename))%>%
    filter(temp_var_ID == varID_of_duplicate)%>%
    select(-temp_var_ID)
  
  merged_as_text<-duplicated_rows %>% 
    select(all_of(collapse_as_text)) %>%
    summarise_all(paste ,collapse=",")
  
  merged_as_numeric_val<-duplicated_rows %>% 
    select(all_of(scoreVals_to_combine)) %>%
    summarise_all(mean)
  
  constant_vals<-duplicated_rows %>% 
    select(all_of(constant_cols)) %>% distinct()

  deduplicated<-cbind(
        constant_vals,
        merged_as_text,
        merged_as_numeric_val
        )%>%
    relocate(any_of(cols_vars))
  if (nrow(deduplicated)>1){
    print("Deduplication not successful!")
    print(constant_vals)
  }
  return(deduplicated)
}

var_ids<-unique(paste(duplicted_vars$ID, duplicted_vars$genename))
#var_ids<-var_ids[1:50]

deduplicated_vars<-lapply(var_ids, deduplicate_vars)

deduplicated_vars_bind<-do.call("rbind", deduplicated_vars)

if (sum(!colnames(deduplicated_vars_bind)==cols_vars)>0) {
  print("Order of columns does not match the requirements!")
}

write_tsv(x=deduplicated_vars_bind,
          file=opt$out_file,
          col_names = FALSE)
