library(readr)
library(tidyverse)
library(data.table)
rm(list=ls())
source("paths_to_folders.R", echo=FALSE)
#820x23
Jolma2013 <- read_delim("../../Data/SELEX-motif-collection/Jolma2013_metadata.csv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE) #22
#593 x 23
Jolma2015 <- read_delim("../../Data/SELEX-motif-collection/Jolma2015_metadata.csv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE) #22

#1 x 23
Morgunova2015 <- read_delim("../../Data/SELEX-motif-collection/Morgunova2015_metadata.csv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE) #22
#10 x 23
Nitta2015 <- read_delim("../../Data/SELEX-motif-collection/Nitta2015_metadata.csv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE) #22

#1161 x 24
Yin2017 <- read_delim("../../Data/SELEX-motif-collection/Yin2017_metadata.csv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE) #23

#1348 x 23
Xie2025 <- vroom::vroom("../../Data/SELEX-motif-collection/Xie2025_metadata.csv", 
                                                     delim = "\t", escape_double = FALSE, 
                                                     trim_ws = TRUE) #22
                         
                         
                         
column_order=c( "ID", "symbol", "Human_Ensemble_ID" ,"clone","family", "Lambert2018_families", 
                "organism", "study","experiment",          
"ligand", "batch","seed", "multinomial","cycle","representative",  "short",               
"type", "comment",  "filename", "IC", "length","consensus","kld_between_revcomp")  

metadata<- rbind( Jolma2013[, column_order], 
                  Jolma2015[, column_order],
                  Nitta2015[, column_order],
                  Morgunova2015[, column_order],
                  Xie2025[,column_order]
                  )

metadata$Methyl.SELEX.Motif.Category=""
metadata=rbind(metadata, Yin2017) #3933

#Families are missing

table(metadata$Lambert2018_families)

#metadata[which(metadata$Lambert2018_families=="Unknown"), 1:7]

metadata[which(metadata$Lambert2018_families=="Unknown"), "Lambert2018_families"]=metadata[which(metadata$Lambert2018_families=="Unknown"), "family"]

tmp=metadata[which(is.na(metadata$Lambert2018_families)), "family"]

#family	Lambert2018_families
# bHLH 	bHLH
# bZIP 	bZIP
# C2H2	C2H2 ZF
# ETS	Ets 
# forkhead	Forkhead
# HMG	HMG
# homeodomain	Homeodomain 
# MEIS	
# nuclearreceptor	Nuclear receptor 
# p53l 	p53
# POU	
# RFX	RFX
# TFAP 	

metadata$family=gsub("C2H2", "C2H2 ZF", metadata$family)
metadata$family=gsub("ETS", "Ets", metadata$family)
metadata$family=gsub("forkhead", "Forkhead", metadata$family)
metadata$family=gsub("homeodomain", "Homeodomain", metadata$family)
metadata$family=gsub("nuclearreceptor", "Nuclear receptor", metadata$family)
metadata$family=gsub("p53l", "p53", metadata$family)

metadata[which(is.na(metadata$Lambert2018_families)), "Lambert2018_families"]=metadata[which(is.na(metadata$Lambert2018_families)), "family"]

#3933 motifs

#Could add better unique identifies
setDT(metadata)

metadata=metadata[,motif_ID := sprintf("ID-%04d", .I)]

# motif_ID column after ID 

metadata <- metadata %>%
  relocate("motif_ID", .after = "ID")


write.table(metadata, file="../../Data/SELEX-motif-collection/metadata.csv", row.names = FALSE, sep="\t")

#Create scdc files from all

#Write .scpd format, all into a single file

append=FALSE

for(i in 1:nrow(metadata)){ #3933
  
  PWM=read.table(paste0( metadata$filename[i]))
  PWM=as.matrix(PWM, dimnames=NULL)
  rownames(PWM)=c("A", "C", "G", "T")
  
  write.table(paste0(">",  metadata$ID[i]),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0(pfms_scpd,"/","all.scpd"))
  append=TRUE
  
  write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0(pfms_scpd,"/","all.scpd"))
  }

line_count <- as.numeric(system(paste("wc -l", paste0(pfms_scpd,"/all.scpd"), intern = TRUE)))
#19665 total This is correct

#Or just concatenate files, something wrong here 

# List all .txt files
#file_list <- list.files(path = "../../Data/PWMs/pfms_scpd", pattern = "\\.scpd$", full.names = TRUE)
#remove all.scpd if it exists

#file_list=file_list[!grepl("all.scpd", file_list)]

# Concatenate and write to output file
#cat(unlist(lapply(file_list, readLines)), sep = "\n", file = paste0(pfms_scpd,"/all_tmp.scpd"))

#line_count <- as.numeric(system(paste("wc -l", paste0(pfms_scpd,"/all_tmp.scpd"), intern = TRUE)))
#13860 total, too less

#save filenames, without "../../

write.table(gsub("../../","",metadata$filename), 
            file="../../Data/SELEX-motif-collection/filenames.csv",
            row.names = FALSE, col.names=FALSE, quote=FALSE)
#save motifnames
write.table(metadata$ID, 
            file="../../Data/SELEX-motif-collection/motifnames.csv",row.names = FALSE, col.names=FALSE, quote=FALSE)


write.table(metadata$motif_ID, 
            file="../../Data/SELEX-motif-collection/motif_IDs.csv",row.names = FALSE, col.names=FALSE, quote=FALSE)

stop()
#Print motif numbers 

metadata %>% filter(study=="Jolma2013"& organism=="Homo_sapiens") %>% nrow() #687
metadata %>% filter(study=="Jolma2013"& organism=="Mus_musculus") %>% nrow() #133

metadata %>% filter(study=="Jolma2015"& experiment=="CAP-SELEX") %>% nrow() #562
metadata %>% filter(study=="Jolma2015"& experiment=="HT-SELEX") %>% nrow() #31
metadata %>% filter(study=="Morgunova2015") %>% nrow() #1
metadata %>% filter(study=="Nitta2015") %>% nrow() #10

metadata %>% filter(study=="Yin2017"& experiment=="HT-SELEX") %>% nrow() #864
metadata %>% filter(study=="Yin2017"& experiment=="Methyl-HT-SELEX") %>% nrow() #297

metadata %>% filter(study=="Xie2025"& experiment=="CAP-SELEX") %>% nrow() #1336
metadata %>% filter(study=="Xie2025"& experiment=="CAP-SELEX") %>% count(type) 
#1: composite  1131
#2:   spacing   205

metadata %>% filter(study=="Xie2025"& experiment=="HT-SELEX") %>% nrow() #12
