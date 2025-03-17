#rm(list=ls())
library(dplyr)
library(readr)
library(stringr)
library(tidyverse)

setwd("/Users/osmalama/projects/TFBS/")

# Reading data, 3933 motifs
df_motif_info <- read_delim('Data/SELEX-motif-collection/metadata.tsv', "\t",
                            escape_double = FALSE, col_names = TRUE, trim_ws = TRUE) 

#1232
sstat_DSA <- read_delim("Results/SELEX-Dominating-Set-Analysis/solutions/SELEX_all_version2.2_min_dom_set_list.txt", 
                        "\t", 
                        col_names = "ID", 
                        escape_double = FALSE, 
                        trim_ws = TRUE)

#The old version of this list has now the wrong motif names, correct



correct_names=bind_rows(lapply(strsplit(sstat_DSA$ID, "_"), function(x) as.data.frame(t(x))))
wrong_inds=which(is.na(correct_names$V8) & !(correct_names$V2 %in% c("HT-SELEX", "Methyl-HT-SELEX")))

wrong_names=sstat_DSA$ID[wrong_inds]
correct_names=correct_names[wrong_inds,]

#There are some with the column 7 value as "NA", move the column values to the right by one starting from the second column
#Do not convert the indexes of the new HT-SELEX motifs

na_ind=which(is.na(correct_names$V7))
correct_names=correct_names[-na_ind,]
wrong_names=wrong_names[-na_ind ]
wrong_inds=wrong_inds[-na_ind]
#correct_names$V7[na_ind]=correct_names$V6[na_ind]
#correct_names$V6[na_ind]=correct_names$V5[na_ind]
#correct_names$V5[na_ind]=correct_names$V4[na_ind]
#correct_names$V4[na_ind]=correct_names$V3[na_ind]
#correct_names$V3[na_ind]=correct_names$V2[na_ind]

#correct_names$V2[na_ind]=NA

#Select the number after m in column 6 and the number after c in column 7 

correct_names$V6=as.numeric(str_extract(correct_names$V6, "\\d+"))

correct_names$V7=as.numeric(str_extract(correct_names$V7, "\\d+"))

# paste rows separated by "_"

correct_names=apply(correct_names, 1, function(x) paste(na.omit(x), collapse = "_"))
names(correct_names)=NULL

sstat_DSA$ID[wrong_inds]=correct_names

# Check for 'representative' column and add 'new_representative' column before that
representative_col_index <- which(colnames(df_motif_info) == "representative")

df_motif_info <- df_motif_info %>%
  add_column(new_representative = "NO", .after = representative_col_index)


#table(is.na(match(sstat_DSA$ID, df_motif_info$ID)))
match_ind=match(sstat_DSA$ID, df_motif_info$ID)
na.ind=which(is.na(match_ind))

head(sstat_DSA$ID[na.ind])

which(df_motif_info$ID==sstat_DSA$ID[na.ind[1]])

#head(sstat_DSA$ID)
#head(df_motif_info$ID[match(sstat_DSA$ID, df_motif_info$ID)])

df_motif_info$new_representative[match(sstat_DSA$ID, df_motif_info$ID)] <- "YES"


# Save the updated dataframe
write_delim(df_motif_info, "Data/SELEX-motif-collection/metadata_representative.tsv", "\t", col_names = TRUE)


# Yimeng's motifs analysis
# 298 + 971=1348
num_yimengs_motifs <- sum(df_motif_info$study == "Xie2025")

df_motif_info %>% filter(study=="Xie2025") %>% count(new_representative)

# new_representative     n
#1 NO                   886
#2 YES                  462

df_motif_info %>% filter(study=="Xie2025") %>% count(type)

# 1 composite  1131
# 2 spacing    205
# 3 HT-selex(NA)   12

df_motif_info %>% filter(study=="Xie2025" & type=="composite") %>% count(new_representative)
#new_representative     n
#NO                   784
#YES                  347

df_motif_info %>% filter(study=="Xie2025" & type=="spacing") %>% count(new_representative)
#new_representative     n
#NO                   93
#YES                  112

df_motif_info %>% filter(study=="Xie2025" & is.na(type)) %>% count(new_representative)
#new_representative     n
#NO                   9
#YES                  3


