rm(list=ls())
library(readr)
metadata_union <- read_delim("/projappl/project_2006203/TFBS/PWMs_final_union/metadata_4003_motifs.csv", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)



metadata <- read_delim("/projappl/project_2006203/TFBS/PWMs_final_version2.2/metadata_3993_motifs.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)



earlier_ID=metadata$ID[which(metadata$ID %in% metadata_union$ID)] #3867
#copy these to new results 

write.table(paste0(earlier_ID, "_top.bed"), file="/projappl/project_2006203/TFBS/PWMs_final_version2.2/MOODS_union.txt",
            col.names=FALSE, row.names = FALSE, quote=FALSE)




new_ID=metadata$ID[which(!(metadata$ID %in% metadata_union$ID))] #66
new_filename=gsub("../../", "../", metadata$filename[which(!(metadata$ID %in% metadata_union$ID))]) #66

write.table(new_ID, file="/projappl/project_2006203/TFBS/PWMs_final_version2.2/new_motifnames.csv",
            col.names=FALSE, row.names = FALSE, quote=FALSE)


write.table(new_filename, file="/projappl/project_2006203/TFBS/PWMs_final_version2.2/new_filenames.csv",
            col.names=FALSE, row.names = FALSE, quote=FALSE)

