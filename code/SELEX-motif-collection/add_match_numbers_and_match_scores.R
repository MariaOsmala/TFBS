
library("readr")
library("tidyverse")
rm(list=ls())
df_motif_info=read_delim("Data/SELEX-motif-collection/metadata_representative.tsv", "\t")


#read match numbers and thresholds 

match_numbers=read_delim("PWMs_final_version2.2/metadata_representatives_match_numbers_thresholds_MeanPhyloP_threshold.tsv", "\t")

correct_names=bind_rows(lapply(strsplit(match_numbers$ID, "_"), function(x) as.data.frame(t(x))))
wrong_inds=which(is.na(correct_names$V8) & !(correct_names$V2 %in% c("HT-SELEX", "Methyl-HT-SELEX", "CAP-SELEX")))

right_inds=which(!is.na(correct_names$V8) | (correct_names$V2 %in% c("HT-SELEX", "Methyl-HT-SELEX", "CAP-SELEX")) )

match_numbers$old_ID=match_numbers$ID
correct_names=correct_names[wrong_inds,]

#There are some with the column 7 value as "NA", move the column values to the right by one starting from the second column


na_ind=which(is.na(correct_names$V7))
#correct_names=correct_names[-na_ind,]

#wrong_inds=wrong_inds[-na_ind]
correct_names$V7[na_ind]=correct_names$V6[na_ind]
correct_names$V6[na_ind]=correct_names$V5[na_ind]
correct_names$V5[na_ind]=correct_names$V4[na_ind]
correct_names$V4[na_ind]=correct_names$V3[na_ind]
correct_names$V3[na_ind]=correct_names$V2[na_ind]

correct_names$V2[na_ind]=NA

#Select the number after m in column 6 and the number after c in column 7 

correct_names$V6=as.numeric(str_extract(correct_names$V6, "\\d+"))

correct_names$V7=as.numeric(str_extract(correct_names$V7, "\\d+"))
correct_names$V8=NULL
correct_names$V9=NULL
# paste rows separated by "_"

correct_names=apply(correct_names, 1, function(x) paste(trimws(na.omit(x)), collapse = "_"))
names(correct_names)=NULL

match_numbers$ID[wrong_inds]=correct_names



match_ind=match(df_motif_info$ID, match_numbers$ID)
na_ind=which(is.na(match_ind))

head(df_motif_info$ID)

head(match_numbers$ID[match_ind])

df_motif_info$match_number <- match_numbers$match_numbers[match_ind]
df_motif_info$threshold <- match_numbers$threshold_double[match_ind]
df_motif_info$phyloP_threshold <- match_numbers$phyloP_threshold[match_ind]

df_motif_info$old_ID=match_numbers$old_ID[match_ind]

#move old_names column after ID 

df_motif_info <- df_motif_info %>% select(ID, old_ID, everything())
#df_motif_info%>% filter(ID==old_ID) %>% nrow()
#2585
#df_motif_info%>% filter(ID!=old_ID) %>% nrow()
#1348

#df_motif_info%>% filter(ID==old_ID) %>% select(study) %>% table()
#df_motif_info%>% filter(ID!=old_ID) %>% select(study) %>% table()

write_delim(df_motif_info, "Data/SELEX-motif-collection/metadata_representative_matchNumber_threshold_phyloP.tsv", "\t", col_names = TRUE)
