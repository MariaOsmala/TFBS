

#############Read excel tables################

library(tidyverse)
library(universalmotif)
rm(list=ls())
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



path="../../PWMs/fromYimeng/PWM_bHLH-Homeobox/"

dir("../../PWMs/fromYimeng/PWM_bHLH-Homeobox/")





  
t="20230420"

#for(file in dir(paste0(path, t, "/"))){
  
#  file.rename(paste0(paste0(path, t, "/"), file), paste0(paste0(path, t, "/"), gsub("-", "_", file ) ) )
  
#} 

for(file in dir(paste0(path, t, "/"))){

    
    
    file.rename(paste0(paste0(path, t, "/"), file), paste0(paste0(path, t, "/"),  gsub("NKX6_1", "NKX6-1", file ) ) )
  
} 



datas=tibble()
data=as_tibble(do.call(rbind,strsplit(dir(paste0(path, t, "/")),"_")))
data$V8=gsub(".pfm", "",as.character(data$V8))
data$symbol=paste(data$V1, data$V2,sep="_")
data=data[,-c(1,2)]



names(data)=c("ligand", "batch", "seed", "multinomial", "cycle", "short", "symbol")
data=data %>% relocate(symbol, .before=ligand)
data$study="fromYimeng_bHLH-Homeobox"
data$type=""
data$organism="Homo_sapiens"
data$clone=NA
data$family=NA
data$experiment="CAP-SELEX"
data$representative=NA
data$comment=NA



data$filename=paste0("PWMs/fromYimeng/PWM_bHLH-Homeobox/", t,"/", dir(paste0(path, t, "/")))
datas=rbind(datas, data)


datas=datas[, c("symbol",	"clone","family",	"organism","study","experiment",
                "ligand",	"batch", "seed",	"multinomial",
                "cycle","representative", "short","type","comment","filename")]

#datas$multinomial=gsub("m", "",as.character(datas$multinomial))
#datas$cycle=gsub("c", "",as.character(datas$cycle))

write.table(datas, file="../../PWMs/fromYimeng/PWM_bHLH-Homeobox_metadata.csv", row.names = FALSE,sep="\t")
saveRDS(datas, file="data/fromYimeng_PWM_bHLH-Homeobox.Rds")

dir.create(paste0("../../PWMs/fromYimeng/pwms_space/PWM_bHLH-Homeobox/"), recursive=TRUE)

dir.create(paste0("../../PWMs/fromYimeng/transfac/PWM_bHLH-Homeobox/"), recursive=TRUE)


#write pwms as .cspd

#write pwms as space separated

append=FALSE

for(m in 1:nrow(datas)){
  f=datas$filename[m]
  PWM=read.table(paste0("../../", f), header=FALSE)
  
  #paste0( strsplit( strsplit(f, "/")[[1]][5], ".pfm")[[1]], "_",datas$type[m])
  write.table(PWM,row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/fromYimeng/pwms_space/",paste0(strsplit(f, "/")[[1]][-c(1,2,4)],collapse="/")) ,sep=" ")
  
  file=paste0("../../PWMs/fromYimeng/pwms_space/",paste0(strsplit(f, "/")[[1]][-c(1,2,4)],collapse="/"))
  
  motif=universalmotif::read_matrix(file=file, sep=" ", header=FALSE)
  motif@name=paste0( strsplit( strsplit(f, "/")[[1]][5], ".pfm")[[1]], "_",datas$type[m])

  transfac=file=paste0("../../PWMs/fromYimeng/transfac/",paste0(strsplit(f, "/")[[1]][-c(1,2,4)],collapse="/"))
  write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
  
  
  #rownames(PWM)=c("A", "C", "G", "T")
  #write.table(paste0(">",  paste0( strsplit( strsplit(f, "/")[[1]][5], ".pfm")[[1]], "_",datas$type[m])),   
  #            append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
  #            file=paste0("../../PWMs/fromYimeng/Homo_sapiens_all", ".scpd"))
  #append=TRUE
  
  #write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
  #            file=paste0("../../PWMs/fromYimeng/Homo_sapiens_all", ".scpd"))
  
}

#display_table_shape(all_table)

