

#############Read excel tables################

library(tidyverse)
library(readxl)
rm(list=ls())

#This table contains background subtracted PWMs for all factors analyzed.  First line for each factor contains information as follows:	
# Symbol:	Fly (Blue colored); Fly gene symbols based on Flybase.org
# Symbol:	Human( Green colored); HGNC symbol (gene name)
# Barcode:	Used barcode
# Batch:	Experimental batch
# Seed:	Sequence used to identify subsequences to be included in the model
# Mul:	Multinomial; Allowed mismatch count to identify subsequences from the seed
# Cyc:	Cycle; SELEX cycle used to construct the model; Previous cycle was used as background unless indicated as "b". Ex) b0 -> Used cycle0 dataset as a background
# "u" means that the model is extracted based on unique reads
# Note	Human NR1I2 and NR1I3 models are not included in the dendrogram analysis (Fig. 2) and Motif network analysis (Fig. S4).

source("code/read_excel_tables.R")



all_table <- read_excel("../../PWMs/Nitta2015/elife-04837-supp1-v1.xlsx", 
                        sheet="E. PWMs", col_names=FALSE, skip=2094)

 
PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(seq(1,nrow(PWMs),5), function(i) PWMs[i,]))

PWMs_metadata=cbind(PWMs_metadata[,1], do.call(rbind, sapply(PWMs_metadata[,2], strsplit, split="_")) )
# "symbol",	"clone","family", "organism",	"study","experiment",
# "ligand",	"batch", "seed",	"multinomial",
# "cycle","representative", "short", "type","comment"

colnames(PWMs_metadata)=c("symbol",	"ligand",	"batch",
                           "seed",	"multinomial",
                          "cycle")

PWMs_metadata$clone=NA
PWMs_metadata$family=NA
PWMs_metadata$organism="Homo_sapiens"
PWMs_metadata$study="Nitta2015"
PWMs_metadata$experiment="HT-SELEX"   
PWMs_metadata$representative=NA
PWMs_metadata$short=NA
PWMs_metadata$type=NA
PWMs_metadata$comment=NA
PWMs_metadata$filename=NA

PWMs_metadata$multinomial=gsub("m", "",as.character(PWMs_metadata$multinomial))
PWMs_metadata$cycle=gsub("c", "",as.character(PWMs_metadata$cycle))

PWMs_metadata=select(PWMs_metadata,c("symbol",	"clone","family", "organism",	"study","experiment",
                                     "ligand",	"batch", "seed",	"multinomial",
                                     "cycle","representative", "short", "type","comment", "filename")
)




PWMs_metadata %>%
  count(symbol) 



PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

dir.create(paste0("../../PWMs/Nitta2015/pwms/","Homo_sapiens"), recursive=TRUE)
dir.create(paste0("../../PWMs/Nitta2015/pwms_space/","Homo_sapiens"), recursive=TRUE)

append=FALSE

for(m in 1:length(PWMs_list)){
  
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Nitta2015/pwms/","Homo_sapiens","/", 
                          paste0(PWMs_metadata[m,
                          -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "short", "type", "representative","filename"))], collapse="_"),
                          ".pfm"), sep="\t")
  
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Nitta2015/pwms_space/","Homo_sapiens","/", 
                          paste0(PWMs_metadata[m,
                                               -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "short", "type", "representative","filename"))], collapse="_"),
                          ".pfm"), sep=" ")
  
  PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
  rownames(PWM)=c("A", "C", "G", "T")
  write.table(paste0(">",  paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone", "family", "organism", "study","comment", "short", "type", "filename"))], collapse="_")),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Nitta2015/",PWMs_metadata[m,"organism"],"_all", ".scpd"))
  append=TRUE
  
  write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Nitta2015/",PWMs_metadata[m,"organism"],"_all", ".scpd"))
  
  PWMs_metadata$filename[m]=paste0("PWMs/Nitta2015/pwms/Homo_sapiens/", paste0(PWMs_metadata[m,
                                                                                           -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "short", "type", "representative","filename"))], collapse="_"),".pfm")
  
  
}


write.table(PWMs_metadata, file="../../PWMs/Nitta2015/metadata.csv", row.names = FALSE, sep="\t")
saveRDS(PWMs_metadata, file="data/Nitta2015.Rds")


#display_table_shape(all_table)

