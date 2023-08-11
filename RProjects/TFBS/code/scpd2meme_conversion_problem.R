

library(readr)
metadata <- read_delim("/projappl/project_2006203/TFBS/PWMs_final/metadata.csv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)



motifs_in_meme=read_delim("/projappl/project_2006203/TFBS/code/motifs_in_meme.txt", delim=",", col_names = FALSE)


metadata[which(!(metadata$ID %in% motifs_in_meme$X1)),1:7]
