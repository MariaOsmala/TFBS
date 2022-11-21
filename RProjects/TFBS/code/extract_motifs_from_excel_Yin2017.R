

#############Read excel tables################

library(tidyverse)
library(readxl)
library(universalmotif)
rm(list=ls())

source("code/read_excel_tables.R")

# This table contains background subtracted PWMs for all of the TF pairs and individual TFs analyzed.  First line for each factor contains information as follows:	
# Base:	A, C, G or T
# symbol(s):	HGNC symbol (gene names)
# clone type:	Type of construct used: extended DNA binding domain (eDBD) or full-length (FL)
# families:	TF structural family
# experiment type:	Type of experiment: HT-SELEX or methyl-HT-SELEX
# ligand sequence:	Name of primer used to generate selection ligand; number before N indicates length of random sequence
# batch:	batch of HT-SELEX and Methyl-HT-SELEX experiment
# seed:	Sequence used to identify subsequences to be included in the model
# multinomial:	Hamming distance to seed used to identify subsequences to be included in the model
# cycle:	SELEX cycle used to construct the model; previous cycle was used as background unless indicated otherwise by a letter "b" followed by the used background cycle, "u" letter indicates that only the first instance of identical sequencing reads were used in the construction of the PWM
# Matrix is one of the representative PWMs:	value yes or no
# comment:	comment field


all_table <- read_excel("../../PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx", 
                        sheet="S2 PWM models", col_names=FALSE, skip=22, n_max=8996-22, col_types=c(rep("text",12), rep("numeric",19)))

 
PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(seq(1,nrow(PWMs),5), function(i) PWMs[i,-1]))



colnames(PWMs_metadata)=c("symbol",	"clone","family",	"experiment",
                          "ligand",	"batch", "seed",	"multinomial",
                          "cycle","representative", "comment")
PWMs_metadata$organism="Homo_sapiens"
PWMs_metadata$study="Yin2017"
PWMs_metadata$short=NA
PWMs_metadata$type=NA


#symbol     clone family organism     study      experiment ligand        batch  seed               multinomial cycle representative short type              comment filename         

PWMs_metadata=select(PWMs_metadata,c("symbol",	"clone","family", "organism",	"study","experiment",
                                     "ligand",	"batch", "seed",	"multinomial",
                                     "cycle","representative", "short", "type","comment")
)


#remove / from family names
PWMs_metadata=PWMs_metadata %>% 
  mutate(family = gsub("/|\\.\\d+[A-Za-z]+", "_aka_", family))

PWMs_metadata %>%
  count(symbol) 

PWMs_metadata$filename=NA


PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

dir.create(paste0("../../PWMs/Yin2017/pwms/","Homo_sapiens"), recursive=TRUE)
dir.create(paste0("../../PWMs/Yin2017/pwms_space/","Homo_sapiens"), recursive=TRUE)
dir.create(paste0("../../PWMs/Yin2017/transfac/","Homo_sapiens"), recursive=TRUE)
append=FALSE

for(m in 1:length(PWMs_list)){
  
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Yin2017/pwms/","Homo_sapiens","/", 
                          paste0(PWMs_metadata[m,
                          -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "short", "type", "filename"))], collapse="_"),
                          ".pfm"), sep="\t")
  
  file=paste0("../../PWMs/Yin2017/pwms/","Homo_sapiens","/", 
              paste0(PWMs_metadata[m,
                                   -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "short", "type", "filename"))], collapse="_"),
              ".pfm")
  
  if ( ncol(PWMs_list[[m]][,-1])!=0 ){
    motif=universalmotif::read_matrix(file=file, sep="\t", header=FALSE)
  
    motif@name=paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone", "family", "organism", "study","comment", "short", "type", "filename"))], collapse="_")
  
    transfac=paste0("../../PWMs/Yin2017/transfac/",PWMs_metadata[m,"organism"],"/", 
                  paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone","family","comment", "study","organism","short", "type","filename"))], collapse="_"),".pfm")
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
  }
  
  
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Yin2017/pwms_space/","Homo_sapiens","/", 
                          paste0(PWMs_metadata[m,
                                               -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "short", "type", "filename"))], collapse="_"),
                          ".pfm"), sep=" ")
  
  PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
  rownames(PWM)=c("A", "C", "G", "T")
  write.table(paste0(">",  paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone", "family", "organism", "study","comment", "short", "type", "filename"))], collapse="_")),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Yin2017/",PWMs_metadata[m,"organism"],"_all", ".scpd"))
  append=TRUE
  
  write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Yin2017/",PWMs_metadata[m,"organism"],"_all", ".scpd"))
  
  PWMs_metadata$filename[m]=paste0("PWMs/Yin2017/pwms/Homo_sapiens/", paste0(PWMs_metadata[m,
                                              -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "short", "type", "filename"))], collapse="_"),".pfm")
  
}

write.table(PWMs_metadata, file="../../PWMs/Yin2017/metadata.csv", row.names = FALSE,sep="\t")
saveRDS(PWMs_metadata, file="data/Yin2017.Rds")




#display_table_shape(all_table)

