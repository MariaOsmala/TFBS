library(readr)

rm(list=ls())
#Lambert families missing
Jolma2013 <- read_delim("~/projects/TFBS/PWMs_final/Jolma2013/metadata.csv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE) #21

Jolma2015 <- read_delim("~/projects/TFBS/PWMs_final/Jolma2015/metadata.csv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE) #22


Morgunova2015 <- read_delim("~/projects/TFBS/PWMs_final/Morgunova2015/metadata.csv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE) #22

colnames(Morgunova2015)[5]="Lambert2018_families"


Nitta2015 <- read_delim("~/projects/TFBS/PWMs_final/Nitta2015/metadata.csv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE) #22

Yin2017 <- read_delim("~/projects/TFBS/PWMs_final/Yin2017/metadata.csv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE) #23

# fromYimeng <- read_delim("~/projects/TFBS/PWMs_final/fromYimeng/metadata.csv", 
#                             delim = "\t", escape_double = FALSE, 
#                             trim_ws = TRUE) #22

fromYimeng_version2 <- read_delim("~/projects/TFBS/PWMs_final_version2/fromYimeng/metadata.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE) #22



column_order=c( "ID", "symbol", "clone","family", "Lambert2018_families", "organism", "study","experiment",          
"ligand", "batch","seed", "multinomial","cycle","representative",  "short",               
"type", "comment",  "filename", "IC", "IC_universal", "length","consensus")  

metadata<- rbind( Jolma2013[, column_order], 
                  Jolma2015[, column_order],
                  Nitta2015[, column_order],
                  Morgunova2015[, column_order],
                  fromYimeng_version2[,column_order]
                  )

metadata$Methyl.SELEX.Motif.Category=""
metadata=rbind(metadata, Yin2017) #3294

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

write.table(metadata, file="../../PWMs_final_version2/metadata.csv", row.names = FALSE, sep="\t")

#Create scdc files from all

#Write .scpd format, all into a single file

append=FALSE

for(i in 1:nrow(metadata)){ #3854
  
  PWM=read.table(paste0( metadata$filename[i]))
  PWM=as.matrix(PWM, dimnames=NULL)
  rownames(PWM)=c("A", "C", "G", "T")
  
  write.table(paste0(">",  metadata$ID[i]),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final_version2/all", ".scpd"))
  append=TRUE
  
  write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final_version2/all", ".scpd"))
  }


#save filenames, without "../../

write.table(gsub("../../","",metadata$filename), file="../../PWMs_final_version2/filenames.csv",row.names = FALSE, col.names=FALSE, quote=FALSE)

write.table(metadata$ID, file="../../PWMs_final_version2/motifnames.csv",row.names = FALSE, col.names=FALSE, quote=FALSE)

#save motifnames

