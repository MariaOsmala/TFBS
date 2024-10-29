library(readr)

library(tidyverse)

library("TFBSTools")
library("universalmotif")
rm(list=ls())

setwd("~/projects/TFBS/RProjects/TFBS/")

metadata=read.csv(file="../../PWMs_final_version2.2/metadata_representatives_match_numbers_thresholds.tsv", sep="\t")

metadata=metadata %>% filter(new_representative=="YES")

i=1
pcm=as.matrix( read.table(paste0( metadata$filename[i]), header=FALSE) )
dimnames(pcm)=list(c("A", "C", "G", "T"))


pcm_class <- TFBSTools::PFMatrix(name=metadata$ID[i], ID=metadata$symbol[i], 
                                 
                                 strand="+", 
                                 bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                 profileMatrix=pcm
)

pwm_class <- TFBSTools::toPWM(pcm_class, type="prob",
                              pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

PPMp <- function(C, p){
  (C + p / nrow(C)) / matrix( (colSums(C) + p),nrow=4, ncol=ncol(C),byrow=TRUE)
  
}

PWM=PPMp(pcm_class@profileMatrix,0.01)

#Position specific scoring matrix
S <- function(C, B) log2(C / B)

PPSM=S(PWM, 0.25)

#Final score for a sequence

score <- function(S, seq) {
  
  sum(S[cbind(match(strsplit(as.character(seq), "")[[1]], rownames(S)), 1:ncol(S))])
}

score(PPSM, "AACGCTAACTAAAATTGGCACAT") # 19.81166

#scan_sequences is from universalmotif
#scan_sequences(PWM,Biostrings::DNAString("AACGCTAACTAAAATTGGCACAT"), threshold.type="logodds")$score #19.8
#scan_sequences(pcm_class,Biostrings::DNAString("AACGCTAACTAAAATTGGCACAT"), threshold.type="logodds")$score #19.801
#scan_sequences(pwm_class@profileMatrix,Biostrings::DNAString("AACGCTAACTAAAATTGGCACAT"), threshold.type="logodds")$score #19.803
#scan_sequences(pwm_class,Biostrings::DNAString("AACGCTAACTAAAATTGGCACAT"), threshold.type="logodds")$score #12.172




#write_homer(PWM, "test.homer", overwrite = TRUE, append = FALSE, 
#            threshold=metadata$threshold_double[i],
#            threshold.type="logodds")


append=FALSE

for(i in 1:nrow(metadata)){ #3854
  
  pcm=as.matrix( read.table(paste0( metadata$filename[i]), header=FALSE) )
  dimnames(pcm)=list(c("A", "C", "G", "T"))
  
  
  pcm_class <- TFBSTools::PFMatrix(name=metadata$ID[i], ID=metadata$symbol[i], 
                                   
                                   strand="+", 
                                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                   profileMatrix=pcm
  )
  
  pwm_class <- TFBSTools::toPWM(pcm_class, type="prob",
                                pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  
  write.table(paste0(">",  metadata$consensus[i], "\t", metadata$ID[i],"\t", metadata$threshold_double[i] ),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final_version2.2/representatives", ".motifs") )
  append=TRUE
  
  write.table(t(pwm_class@profileMatrix),  append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final_version2.2/representatives", ".motifs"), sep="\t")
  
  
  
  
  
}





# 
# 
# 
# Arguments
# motifs See convert_motifs() for acceptable formats.
# file character(1) File name.
# overwrite logical(1) Overwrite existing file.
# append logical(1) Add to an existing file.
# 
# ?write_homer




# # Open a file connection
# file_conn <- file(file_name, "w")
# 
# # Write the identifier and description in JASPAR format
# writeLines(paste0(">", identifier, " ", description), file_conn)
# 
# # Write each nucleotide and its counts
# for (base in rownames(matrix)) {
#   counts <- matrix[base, ]
#   counts_str <- paste(counts, collapse = "  ")
#   writeLines(paste(base, "[", counts_str, "]"), file_conn)
# }
# 
# # Close the file connection
# close(file_conn)
# 
# append=FALSE
# 
# for(i in 1:nrow(metadata)){ #3854
#   
#   PWM=read.table(paste0( metadata$filename[i]))
#   PWM=as.matrix(PWM, dimnames=NULL)
#   rownames(PWM)=c("A", "C", "G", "T")
#   
#   write.table(paste0(">",  metadata$ID[i], "\t", metadata$symbol[i]),   
#               append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
#               file=paste0("../../PWMs_final_version2.2/representatives", ".jaspar"))
#   append=TRUE
#   
#   # Write each nucleotide and its counts
#   for (base in rownames(PWM)) {
#     counts <- PWM[base, ]
#     counts_str <- paste(counts, collapse = " ")
#     
#     write.table(paste0(base," ", "[ ", counts_str, " ]") ,append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
#                 file=paste0("../../PWMs_final_version2.2/representatives", ".jaspar"))
#     
#   }
#   
#   
#  
#   }
# 

