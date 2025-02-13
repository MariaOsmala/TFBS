

#############Read excel tables################

library(tidyverse)
library(readxl)
library(universalmotif)
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

pwms="../../PWMs_final/Nitta2015/pwms/Homo_sapiens"
pwms_space="../../PWMs_final/Nitta2015/pwms_space/Homo_sapiens"
pwms_transfac="../../PWMs_final/Nitta2015/transfac/Homo_sapiens"

dir.create(pwms, recursive=TRUE)
dir.create(pwms_space, recursive=TRUE)
dir.create(pwms_transfac, recursive=TRUE)

append=FALSE

#Remove these from metadata
remove=c()

for(m in 1:length(PWMs_list)){
  
  #Empty matrix
  if( ncol(PWMs_list[[m ]]) ==1 ) {
    remove=c(remove, m)
    print(m)
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWMs_list[[ m ]][,-1]))) ==0 ){
    remove=c(remove, m)    
    print(m)
  }else{
    
    filename=paste0(pwms,"/", 
                    paste0(PWMs_metadata[m,
                                         -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", 
                                                                              "representative", "short", "type", "filename","ID", "IC", "IC_universal","length", "consensus"))], collapse="_"),
                    ".pfm")
    ID=paste0(PWMs_metadata[m,
                            -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative", "short", "type", "filename",
                                                                 "ID", "IC", "IC_universal","length", "consensus" ))], collapse="_")
    #Filename and ID should be unique
    if(filename %in% PWMs_metadata$filename){
      print("Warning! Non-unique filenames and IDs")
    }
    
    PWMs_metadata$filename[m]=filename
    PWMs_metadata$ID[m]=ID
  
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=filename, sep="\t")
  
    pcm=as.matrix( read.csv(filename, header=FALSE, sep="\t") )
    colnames(pcm)=NULL
    rownames(pcm)=c("A", "C", "G", "T")
    
    pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                profileMatrix=pcm
    )
    
    #pwm <- TFBSTools::toPWM(pfm, type="log2probratio", pseudocounts=0.01)
    icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
    PWMs_metadata[m, "IC"]=sum(rowSums(icm)) #6.77892 # universalmotif computes this wrong?
    
    #motif length
    PWMs_metadata[m, "length"]=length(pfm)
    motif=universalmotif::read_matrix(file=filename, sep="\t", header=FALSE)
    motif@name=ID
    PWMs_metadata[m, "IC_universal"]=motif@icscore
    PWMs_metadata[m, "consensus"]=motif@consensus
  
    transfac=paste0(pwms_transfac, "/", ID,".pfm")
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
  
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0(pwms_space,"/", ID, ".pfm"), sep=" ")
  
  PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
  rownames(PWM)=c("A", "C", "G", "T")
  write.table(paste0(">",  ID),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final/Nitta2015/all", ".scpd"))
  append=TRUE
  
  write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final/Nitta2015/all", ".scpd"))
  
  }
}

TF_family_mappings <- read_csv("~/projects/motif-clustering-Viestra-private/HumanTFs-ccbr-TableS1.csv")


# Merge two data frames
PWMs_metadata <- PWMs_metadata %>%
  left_join(select(TF_family_mappings, "symbol", "Gene Information/ID","DBD" ), by = 'symbol')

PWMs_metadata <- add_column(PWMs_metadata, Lambert2018_families = PWMs_metadata$DBD, .before = 4)
PWMs_metadata$`Gene Information/ID`=NULL
PWMs_metadata$DBD=NULL



write.table(PWMs_metadata, file="../../PWMs_final/Nitta2015/metadata.csv", row.names = FALSE, sep="\t")
saveRDS(PWMs_metadata, file="Rdata/Nitta2015.Rds")


#display_table_shape(all_table)

