

#############Read excel tables################

library(tidyverse)
library(universalmotif)

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)



rm(list=ls())

types=c("Batch_spacing_final_v2", "BatchOverlap_final_v2") #spacing composite

path="../../PWMs_final_version2/fromYimeng/"


############## Remove dots from the filename ################
# Set the target directory
# directory <- paste0(path, types[1]) # Change this to your directory path
# directory <- paste0(path, types[2]) # Change this to your directory path
# # List all files in the directory
# files <- list.files(directory, full.names = TRUE)
# 
# # Loop through each file to rename it
# for (file in files) {
#   #file=files[1]
#   # Extract the directory, filename, and extension
#   dir_name <- dirname(file)
#   base_name <- basename(file)
#   file_ext <- tools::file_ext(base_name)
#   name_without_ext <- tools::file_path_sans_ext(base_name)
#   
#   # Replace dots in the filename, excluding the extension
#   new_name_without_ext <- gsub("\\.", "_", name_without_ext)
#   
#   # Construct the new filename with the original extension
#   new_filename <- ifelse(nchar(file_ext) > 0, 
#                          paste0(new_name_without_ext, ".", file_ext), 
#                          new_name_without_ext)
#   
#   # Construct the full path for the new filename
#   new_file_path <- file.path(dir_name, new_filename)
#   
#   # Rename the file
#   file.rename(file, new_file_path)
# }
# 



datas=tibble()
  
for(t in types){
  #t=types[1]
  data=as_tibble(do.call(rbind,strsplit(dir(paste0(path, t, "/")),split="_")))
  data$V8=gsub(".pfm", "", data$V8)
  data$symbol=paste(data$V1, data$V2,sep="_")
  data=data[,-c(1,2)]
  names(data)=c("ligand", "batch", "seed", "multinomial", "cycle", "short", "symbol")
  data=data %>% relocate(symbol, .before=ligand)
  data$study="fromYimeng_version2"
  data$type=t
  data$organism="Homo_sapiens"
  data$clone=NA
  data$family=NA
  data$experiment="CAP-SELEX"
  data$representative=NA
  data$comment=NA
  data$filename=paste0("PWMs_final_version2/fromYimeng/", t,"/", dir(paste0(path, t, "/")))
  datas=rbind(datas, data)
}

datas=datas[, c("symbol",	"clone","family",	"organism","study","experiment",
                "ligand",	"batch", "seed",	"multinomial",
                "cycle","representative", "short","type","comment","filename")]

pwms=paste0("../../PWMs_final_version2/fromYimeng/pwms/Homo_sapiens/")
pwms_space=paste0("../../PWMs_final_version2/fromYimeng/pwms_space/Homo_sapiens/")
transfac_path=paste0("../../PWMs_final_version2/fromYimeng/transfac/Homo_sapiens/")
dir.create(pwms, recursive=TRUE)
dir.create(pwms_space, recursive=TRUE)
dir.create(transfac_path, recursive=TRUE)

#write pwms as .cspd
#write pwms as space separated

append=FALSE

datas$ID=""
datas$filename_final=""
datas$IC=NA
datas$length=NA
datas$IC_universal=NA
datas$consensus=""

#Remove these from metadata, empty TFs
remove=c()

table(datas$type)
# Batch_spacing_final_v2  BatchOverlap_final_v2 
# 298                    971 

length(unique(datas %>%
  rowwise() %>%
  mutate(concatenated = paste(symbol, ligand, batch, seed, multinomial, cycle, sep = "_")) %>% select(concatenated) %>%pull(concatenated) )) #These are all unique

for(m in 1:nrow(datas)){
  #m=1#
  f=datas$filename[m]
  PWM=read.table(paste0("../../", f))
  #Empty matrix
  if( ncol(PWM) ==1 ) {
    remove=c(remove, m)
    print(f)
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWM[,-1]))) ==0 ){
    remove=c(remove, m)    
    print(f)
  }else{
  
    filename=paste0(strsplit(f, "/")[[1]][-c(1,2,3)],collapse="/")
    #add extension
    filename_final=paste0(gsub(".pfm", "", filename) ,".pfm")
    
    
    if(paste0(pwms,filename_final) %in% datas$filename_final){
      print("Warning! Non-unique filenames and IDs")
      print(m)
    }
    
  
    ID=gsub(".pfm", "",filename_final)
    datas$ID[m]=ID
    
    #write if to pwms, tab separated
    write.table(PWM,row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(pwms,filename_final) ,sep="\t")
    
    
    
    #paste0( strsplit( strsplit(f, "/")[[1]][5], ".pfm")[[1]], "_",datas$type[m])
    write.table(PWM,row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(pwms_space,filename_final) ,sep=" ")
    
    
    pcm=as.matrix( read.csv(paste0(pwms,filename_final), header=FALSE, sep="\t") )
    colnames(pcm)=NULL
    rownames(pcm)=c("A", "C", "G", "T")
    
    pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                profileMatrix=pcm
    )
    
    #pwm <- TFBSTools::toPWM(pfm, type="log2probratio", pseudocounts=0.01)
    icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
   
    datas$IC[m]=sum(rowSums(icm)) #6.77892 # universalmotif computes this wrong?
    
    #motif length
    datas$length[m]=length(pfm)
    
    motif=universalmotif::read_matrix(file=paste0(pwms,filename_final), sep="\t", header=FALSE)
    
   
    datas[m, "IC_universal"]=motif@icscore
    datas[m, "consensus"]=motif@consensus
    
    
    #Write transfac format
    
    motif=universalmotif::read_matrix(file=paste0(pwms_space,filename_final), sep=" ", header=FALSE)
    motif@name=ID
  
    transfac_file=paste0(transfac_path,filename_final)
    write_transfac(motif, file=transfac_file, overwrite = TRUE, append = FALSE)
    
    #Write .scpd format, all into a single file
    PWM=as.matrix(PWM, dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
                append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0("../../PWMs_final_version2/fromYimeng/Homo_sapiens_all", ".scpd"))
    append=TRUE
    
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
                file=paste0("../../PWMs_final_version2/fromYimeng/Homo_sapiens_all", ".scpd"))
    
    
    datas$filename_final[m]=paste0(pwms,filename_final)
    
  
  }
}

#There is nothing to be removed
#datas=datas[-remove,]

#Add family mappings

TF_family_mappings <- read_delim("~/projects/motif-clustering-Viestra-private/HumanTFs-ccbr-TableS1.csv", delim=",")

pairs <- as.data.frame(do.call(rbind, strsplit(datas$symbol, split="_")))
colnames(pairs)=c("symbol", "symbol1")
pairs$order=seq(nrow(pairs))

commondf_pairs <- merge(pairs, TF_family_mappings, by = 'symbol', all.x = TRUE)[,c(1,2,3,5)]
commondf_pairs <- commondf_pairs[order(commondf_pairs$order), ]
commondf_pairs$order <- NULL

names(commondf_pairs)=c("symbol0", "symbol", "DBD0")
commondf_pairs$order=seq(nrow(commondf_pairs))
commondf_pairs2 = merge(commondf_pairs, TF_family_mappings, by = 'symbol', all.x = TRUE)[,c(1,2,3,4,6)]
commondf_pairs2 <- commondf_pairs2[order(commondf_pairs2$order), ]
commondf_pairs2$order <- NULL
commondf_pairs2 <- commondf_pairs2[, c(2,1,3,4)]
names(commondf_pairs2)=c("symbol0", "symbol1", "DBD0", "DBD1")

#Are there some TFs for which there is no family info
unique(commondf_pairs2[which(is.na(commondf_pairs2$DBD0)),"symbol0"])
# "BHLHB2" bHLH
# "BHLHB5" bHLH
# "BHLHB8" bHLH
# "gntR"  GntR-type HTH domain 
# "FOXO3A" "Forkhead
# "ZNF238" "C2H2 ZF"

commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB2")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB5")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB8")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="gntR")]="GntR-type HTH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="FOXO3A")]="Forkhead"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="ZNF238")]="C2H2 ZF"

table(commondf_pairs2$DBD0, useNA="always") #There are no NAs
table(commondf_pairs2$DBD1, useNA="always") #There are no NAs


pair_families <- paste(commondf_pairs2$DBD0,commondf_pairs2$DBD1, sep = "_")

datas <- add_column(datas, Lambert2018_families = pair_families, .before = 4)

unique_Lambertfamilies <- unique(datas$Lambert2018_families)

datas$filename=NULL
datas <- rename(datas, filename = filename_final)

write.table(datas, file="../../PWMs_final_version2/fromYimeng/metadata.csv", row.names = FALSE,sep="\t")
saveRDS(datas, file="RData/fromYimeng_version2.Rds")


