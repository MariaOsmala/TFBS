#rm(list=ls())
library(dplyr)
library(readr)
library(stringr)
library(tidyverse)

setwd("/Users/osmalama/projects/TFBS/RProjects/TFBS")

# Reading data, 3933!
df_motif_info <- read_delim('../../PWMs_final_version2.2/metadata_3993_motifs.csv', "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE) 

new_composites=df_motif_info %>% filter(study=="fromYimeng" & type=="composite") #1131

#391
sstat_DSA <- read_delim("../../../SELEX-Dominating-Set-Analysis/solutions/SELEX_new_composites_version2.2_min_dom_set_list.txt", "\t", col_names = "ID", escape_double = FALSE, trim_ws = TRUE)

# Check for 'representative' column and add 'new_representative' column before that
representative_col_index <- which(colnames(new_composites) == "representative")

new_composites <- new_composites %>%
  add_column(new_representative = "NO", .after = representative_col_index)

new_composites$new_representative[match(sstat_DSA$ID, new_composites$ID)] <- "YES"


# Save the updated dataframe
write_delim(new_composites, "../../PWMs_final_version2.2/new_composites_representatives.tsv", "\t", col_names = TRUE)


# Yimeng's motifs analysis
# 298 + 971=1348
num_yimengs_motifs <- sum(df_motif_info$study == "fromYimeng")

new_composites %>% count(new_representative)

# new_representative     n
#1 NO                   740
#2 YES                  391

