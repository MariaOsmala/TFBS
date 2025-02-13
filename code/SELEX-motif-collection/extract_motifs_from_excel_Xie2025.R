
library(tidyverse)
library(universalmotif)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)



rm(list=ls())

#Version2.2

#HT-SELEX
length(dir("/Users/osmalama/projects/TFBS/PWMs_final_version2.2/novel_motif_HT-SELEX/")) #12

path="../../PWMs_final_version2.2/novel_motif_HT-SELEX/"

datas=tibble()
data=as_tibble(do.call(rbind,strsplit(dir(path),split="_")))
data$V7=gsub(".pfm", "", data$V7)
names(data)=c("symbol","ligand", "batch", "seed", "multinomial", "cycle", "short")
data$study="fromYimeng"
data$type=NA
data$organism="Homo_sapiens"
data$clone=NA
data$family=NA
data$experiment="HT-SELEX"
data$representative=NA
data$comment=NA
data$filename=paste0("PWMs_final_version2.2/novel_motif_HT-SELEX/", dir(paste0(path)))
data$ID=gsub(".pfm", "",do.call(rbind, strsplit(data$filename, "/"))[,3])
data=data %>% relocate(ID, .before=symbol)
datas=rbind(datas, data)




datas=datas[, c("ID","symbol",	"clone","family",	"organism","study","experiment",
                "ligand",	"batch", "seed",	"multinomial",
                "cycle","short","representative", "type","comment","filename")]






types=c("Composite_motif_version2.2","Spacing_motif_version2.2")


path="../../PWMs_final_version2.2/fromYimeng/"



for(t in types){
  
  #t=types[1]
  data=as_tibble(do.call(rbind,strsplit(dir(paste0(path, t, "/")),split="_")))
  data$V7=gsub(".pfm", "", data$V7)
  data$symbol=paste(data$V1, data$V2,sep="_")
  data=data[,-c(1,2)]
  names(data)=c("ligand", "batch", "seed", "multinomial", "cycle", "symbol")
  data=data %>% relocate(symbol, .before=ligand)
  data$short=NA
  data$study="fromYimeng"
  data$type=t #composite   spacing   unclear
  data$organism="Homo_sapiens"
  data$clone=NA
  data$family=NA
  data$experiment="CAP-SELEX"
  data$representative=NA
  data$comment=NA
  data$filename=paste0("PWMs_final_version2.2/fromYimeng/", t,"/", dir(paste0(path, t, "/")))
  data$ID=gsub(".pfm", "",do.call(rbind, strsplit(data$filename, "/"))[,4])
  data=data %>% relocate(ID, .before=symbol)
  data=data[, c("ID","symbol",	"clone","family",	"organism","study","experiment",
                  "ligand",	"batch", "seed",	"multinomial",
                  "cycle","short","representative", "type","comment","filename")]
  
  datas=rbind(datas, data)
}







pwms=paste0("../../PWMs_final_version2.2/fromYimeng/pwms/Homo_sapiens/")
pwms_space=paste0("../../PWMs_final_version2.2/fromYimeng/pwms_space/Homo_sapiens/")
transfac_path=paste0("../../PWMs_final_version2.2/fromYimeng/transfac/Homo_sapiens/")
dir.create(pwms, recursive=TRUE)
dir.create(pwms_space, recursive=TRUE)
dir.create(transfac_path, recursive=TRUE)

#write pwms as .cspd
#write pwms as space separated

append=FALSE

datas$filename_final=""
datas$IC=NA
datas$length=NA
datas$IC_universal=NA
datas$consensus=""
#Remove these from metadata, empty TFs
remove=c()


for(m in 1:nrow(datas)){
  print(m)
  f=datas$filename[m]
  PWM=as.matrix(read_tsv(paste0("../../", f), col_names=FALSE))
  #Empty matrix
  if( ncol(PWM) ==1 ) {
    remove=c(remove, m)
    print(paste0("Empty matrix: ", f))
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWM[,-1]))) ==0 ){
    remove=c(remove, m)    
    print(paste0("Zero matrix: ",f))
  }else{
  
    if(datas$experiment[m]=="HT-SELEX"){
      filename=paste0(strsplit(f, "/")[[1]][-c(1,2)],collapse="/")  
    }else{
      #CAP-SELEX
      filename=paste0(strsplit(f, "/")[[1]][-c(1,2,3)],collapse="/")
    }
    
    #add extension
    if(length(grep("_short.pfm", filename))!=0){
      filename_final=paste0(gsub("_short.pfm", "", filename), ".pfm") 
    }else{
      filename_final=filename
    }
    
    
  
    ID=gsub(".pfm", "",filename_final)
    datas$ID[m]=ID
    if(!is.na(datas$type[m])){
      if(datas$type[m]=="Composite_motif_version2.2"){
        datas$type[m]="composite"
      }else if(datas$type[m]=="Spacing_motif_version2.2"){
        datas$type[m]="spacing"
      }
    }
    
    
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
    
    transfac=paste0(transfac_path,filename_final)
    
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
    
    
    motif=universalmotif::read_matrix(file=paste0(pwms_space,filename_final), sep=" ", header=FALSE)
    motif@name=ID
  
    transfac_file=paste0(transfac_path,filename_final)
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
    #Write .scpd format, all into a single file
    PWM=as.matrix(PWM, dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
                append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0("../../PWMs_final_version2.2/fromYimeng/Homo_sapiens_all", ".scpd"))
    append=TRUE
    
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
                file=paste0("../../PWMs_final_version2.2/fromYimeng/Homo_sapiens_all", ".scpd"))
    
    
    datas$filename_final[m]=paste0(pwms,filename_final)
    
  
  }
}



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
#"THHEX"
#"BHLHB2" bHLH
#"BHLHB8" bHLH
#"FOXO3A" Forkhead
# "gntR"  GntR-type HTH domain 
# "ZNF238" "C2H2 ZF"
#"BHLHB5" bHLH

commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB2")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB5")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB8")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="FOXO3A")]="Forkhead"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="gntR")]="GntR-type HTH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="ZNF238")]="C2H2 ZF"




table(commondf_pairs2$DBD0, useNA="always") #There is 1 NA
table(commondf_pairs2$DBD1, useNA="always") #There is 1 NA


pair_families <- paste(commondf_pairs2$DBD0,commondf_pairs2$DBD1, sep = "_")

#only one family for HT-SELEX motifs

pair_families[which(datas$experiment=="HT-SELEX")]=do.call(rbind, strsplit(pair_families[which(datas$experiment=="HT-SELEX")],"_"))[,1]


datas <- add_column(datas, Lambert2018_families = pair_families, .before = 4)
unique_Lambertfamilies <- unique(datas$Lambert2018_families)

datas$filename=NULL
datas <- rename(datas, filename = filename_final)

column_order=c( "ID", "symbol", "clone","family", "Lambert2018_families", "organism", "study","experiment",          
                "ligand", "batch","seed", "multinomial","cycle","representative",  "short",               
                "type", "comment",  "filename", "IC", "IC_universal", "length","consensus")  

datas=datas[,column_order]



write.table(datas, file="../../PWMs_final_version2.2/fromYimeng/metadata.csv", row.names = FALSE,sep="\t")

saveRDS(datas, file="RData/fromYimeng_version2.2.Rds")
datas=readRDS(file="RData/fromYimeng_version2.2.Rds")

table(datas$experiment)
#CAP-SELEX  HT-SELEX 
# 1336        12 

table(datas$type)
# composite   spacing 
# 1131       205 

#Check that the pfms match

datas_union=readRDS(file="RData/fromYimeng_union.Rds")
missing_in_union=c()

for(m in 1:nrow(datas)){
  version2.2_pfm=read.table(datas$filename[m]) #.pfm is missing from the filename?
  
  match_ind=which(datas_union$ID==datas$ID[m])
  if(length(match_ind)==0){
    missing_in_union=c(missing_in_union, m)
    #print(m)
    #print("ID not found in union")
    
  }else{
    
    version_union_pfm=read.table(paste0(datas_union$filename[match_ind]))
    if( is.na(table(version2.2_pfm==version_union_pfm)["FALSE"]) ){
      #Same motifs!
    }else{
      print(m)
      print("pfms of unique IDs are not the same!")
    }
    
  }
}

#Which are missing in union (66)
datas$ID[missing_in_union]

write.table(datas$ID[missing_in_union], file="../../PWMs_final_version2.2/fromYimeng/motifs_missing_in_union.csv", row.names = FALSE,sep="\t")
write.table(gsub("../../","", datas$filename[missing_in_union]), file="../../PWMs_final_version2.2/fromYimeng/filenames_missing_in_union.csv", row.names = FALSE,sep="\t")
