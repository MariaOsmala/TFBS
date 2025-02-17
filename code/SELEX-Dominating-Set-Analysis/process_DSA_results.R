#rm(list=ls())
library(dplyr)
library(readr)
library(stringr)
library(tidyverse)

setwd("/Users/osmalama/projects/TFBS/RProjects/TFBS")

# Reading data, 3933!
df_motif_info <- read_delim('../../PWMs_final_version2.2/metadata_3993_motifs.csv', "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE) 

#1232
sstat_DSA <- read_delim("../../../SELEX-Dominating-Set-Analysis/solutions/SELEX_all_version2.2_min_dom_set_list.txt", "\t", col_names = "ID", escape_double = FALSE, trim_ws = TRUE)

# Check for 'representative' column and add 'new_representative' column before that
representative_col_index <- which(colnames(df_motif_info) == "representative")

df_motif_info <- df_motif_info %>%
  add_column(new_representative = "NO", .after = representative_col_index)



#table(is.na(match(sstat_DSA$ID, df_motif_info$ID)))
#head(sstat_DSA$ID)
#head(df_motif_info$ID[match(sstat_DSA$ID, df_motif_info$ID)])

df_motif_info$new_representative[match(sstat_DSA$ID, df_motif_info$ID)] <- "YES"


# Save the updated dataframe
write_delim(df_motif_info, "../../PWMs_final_version2.2/metadata_representatives.tsv", "\t", col_names = TRUE)


# Yimeng's motifs analysis
# 298 + 971=1348
num_yimengs_motifs <- sum(df_motif_info$study == "fromYimeng")

df_motif_info %>% filter(study=="fromYimeng") %>% count(new_representative)

# new_representative     n
#1 NO                   886
#2 YES                  462

df_motif_info %>% filter(study=="fromYimeng") %>% count(type)

# 1 composite  1131
# 2 spacing    205
# 3 HT-selex   12

df_motif_info %>% filter(study=="fromYimeng" & type=="composite") %>% count(new_representative)
#new_representative     n
#NO                   784
#YES                  347

df_motif_info %>% filter(study=="fromYimeng" & type=="spacing") %>% count(new_representative)
#new_representative     n
#NO                   93
#YES                  112

df_motif_info %>% filter(study=="fromYimeng" & is.na(type)) %>% count(new_representative)
#new_representative     n
#NO                   9
#YES                  3

###### Table for google drive ######

#fromYimeng_metadata <- read_delim("~/projects/TFBS/PWMs_final_version2.2/fromYimeng/metadata.csv", #1348
#     delim = "\t", escape_double = FALSE, 
#     trim_ws = TRUE)

all_motifs_metadata <- read_delim("~/projects/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", #1348
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

fromYimeng_metadata=all_motifs_metadata %>% filter(study=="fromYimeng")



fromYimeng_metadata=fromYimeng_metadata %>% select(c("ID", "experiment", "new_representative", "type"          ))

#rename columns
colnames(fromYimeng_metadata) <- c("ID", "experiment", "representative", "type")

write_delim(fromYimeng_metadata, "../../PWMs_final_version2.2/share_google_drive/representatives.tsv", "\t", col_names = TRUE)



Supplementary_Table_S2_PWMs <- read_csv("~/projects/TFBS/PWMs_final_version2.2/share_google_drive/Supplementary_Tables - Table S2_PWMs.csv", col_names = FALSE)

nrow(Supplementary_Table_S2_PWMs)/5 # 1269

