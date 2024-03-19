
library(tidyverse)
library(universalmotif)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)



rm(list=ls())



# Version 1 motifs --------------------------------------------------------

types=c("pfm_composite_new", "pfm_spacing_new", "20230420")

path="../../PWMs/fromYimeng/pfm_from_newData/"

datas=tibble()
  
for(t in types){
   #t=types[1]
  data=as_tibble(do.call(rbind,strsplit(dir(paste0(path, t, "/")),"_")))
  data$V8=gsub(".pfm", "",as.character(data$V8))
  data$symbol=paste(data$V1, data$V2,sep="_")
  data=data[,-c(1,2)]
  names(data)=c("ligand", "batch", "seed", "multinomial", "cycle", "short", "symbol")
  data=data %>% relocate(symbol, .before=ligand)
  data$study="fromYimeng"
  data$type=t
  data$organism="Homo_sapiens"
  data$clone=NA
  data$family=NA
  data$experiment="CAP-SELEX"
  data$representative=NA
  data$comment=NA
  data$filename=paste0("PWMs/fromYimeng/pfm_from_newData/", t,"/", dir(paste0(path, t, "/")))
  datas=rbind(datas, data)
}

datas=datas[, c("symbol",	"clone","family",	"organism","study","experiment",
                "ligand",	"batch", "seed",	"multinomial",
                "cycle","representative", "short","type","comment","filename")]

#datas$multinomial=gsub("m", "",as.character(datas$multinomial))
#datas$cycle=gsub("c", "",as.character(datas$cycle))

#write.table(datas, file="../../PWMs/fromYimeng/metadata.csv", row.names = FALSE,sep="\t")
#saveRDS(datas, file="data/fromYimeng.Rds")



pwms=paste0("../../PWMs_final_union/fromYimeng/pwms/Homo_sapiens/")
pwms_space=paste0("../../PWMs_final_union/fromYimeng/pwms_space/Homo_sapiens/")
transfac_path=paste0("../../PWMs_final_union/fromYimeng/transfac/Homo_sapiens/")
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

#capselex=readRDS(file="/Users/osmalama/projects/TFBS/RProjects/TFBS/RData/half_site_recognition_in_capselex_motifs.RDS")
#test=capselex[grep("HOXB13_MEIS3_TGCGAT40NGAT_YRIIII_MATAAANNTGTCA_m1_c3b0u", capselex$ID),]


datas$filename[grep("HOXB13_MEIS3_TGCGAT40NGAT_YRIIII_MATAAANNTGTCA_m1_c3b0u", datas$filename)]
#The .pfms are equal, set this to both types: pfm_composite_new", "pfm_spacing_new"
# "PWMs/fromYimeng/pfm_from_newData/pfm_composite_new/HOXB13_MEIS3_TGCGAT40NGAT_YRIIII_MATAAANNTGTCA_m1_c3b0u_short.pfm"
# "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/HOXB13_MEIS3_TGCGAT40NGAT_YRIIII_MATAAANNTGTCA_m1_c3b0u_short.pfm" 

datas$filename[grep("HOXD10_TBX6_TCTACT40NACA_YAAII_NYMRTAAAANAGGTGTNA_m1_c3b0", datas$filename)]
#"PWMs/fromYimeng/pfm_from_newData/pfm_composite_new/HOXD10_TBX6_TCTACT40NACA_YAAII_NYMRTAAAANAGGTGTNA_m1_c3b0_short.pfm"
#"PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/HOXD10_TBX6_TCTACT40NACA_YAAII_NYMRTAAAANAGGTGTNA_m1_c3b0_short.pfm"  

both_types=c("HOXB13_MEIS3_TGCGAT40NGAT_YRIIII_MATAAANNTGTCA_m1_c3b0u", "HOXD10_TBX6_TCTACT40NACA_YAAII_NYMRTAAAANAGGTGTNA_m1_c3b0")

#Remove the other instance of these from the datas

datas=datas[-grep(both_types[1], datas$filename)[1],]
datas=datas[-grep(both_types[2], datas$filename)[1],]

for(m in 1:nrow(datas)){
  f=datas$filename[m]
  PWM=read.table(paste0("../../", f))
  #Empty matrix
  if( ncol(PWM) ==1 ) {
    remove=c(remove, m)
    print(paste0("Empty matrix: ", f))
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWM[,-1]))) ==0 ){
    remove=c(remove, m)    
    print(paste0("Zero matrix: ",f))
  }else{
  
    filename=paste0(strsplit(f, "/")[[1]][-c(1,2,3,4)],collapse="/")
    #add extension
    filename_final=paste0(gsub("_short.pfm", "", filename), ".pfm") 
    
    
    if(paste0(pwms,filename_final) %in% datas$filename_final){
      print("Warning! Non-unique filenames and IDs")
      print(filename_final)
    }
    
  
    ID=gsub(".pfm", "",filename_final)
    datas$ID[m]=ID
    
    if(datas$type[m]=="pfm_composite_new"){
      datas$type[m]="composite"
    }else if(datas$type[m] %in% c("pfm_spacing_new", "20230420")){
      datas$type[m]="spacing"
    }
    
    if(ID %in% both_types){
      datas$type[m]="both"
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
    
    #file=paste0("../../PWMs/fromYimeng/pwms_space/",paste0(strsplit(f, "/")[[1]][-c(1,2,3)],collapse="/"))
    
    motif=universalmotif::read_matrix(file=paste0(pwms_space,filename_final), sep=" ", header=FALSE)
    motif@name=ID
  
    transfac_file=paste0(transfac_path,filename_final)
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
    #Write .scpd format, all into a single file
    PWM=as.matrix(PWM, dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
                append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0("../../PWMs_final_union/fromYimeng/Homo_sapiens_all", ".scpd"))
    append=TRUE
    
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
                file=paste0("../../PWMs_final_union/fromYimeng/Homo_sapiens_all", ".scpd"))
    
    
    datas$filename_final[m]=paste0(pwms,filename_final)
    
  
  }
}

datas=datas[-remove,] #707

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
#"BHLHB2" bHLH
#"BHLHB8" bHLH
#"FOXO3A" Forkhead
#"BHLHB5" bHLH

commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB2")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB5")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB8")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="FOXO3A")]="Forkhead"

table(commondf_pairs2$DBD0, useNA="always") #There are no NAs
table(commondf_pairs2$DBD1, useNA="always") #There are no NAs


pair_families <- paste(commondf_pairs2$DBD0,commondf_pairs2$DBD1, sep = "_")



datas <- add_column(datas, Lambert2018_families = pair_families, .before = 4)
unique_Lambertfamilies <- unique(datas$Lambert2018_families)

datas$filename=NULL
datas <- rename(datas, filename = filename_final)

column_order=c( "ID", "symbol", "clone","family", "Lambert2018_families", "organism", "study","experiment",          
                "ligand", "batch","seed", "multinomial","cycle","representative",  "short",               
                "type", "comment",  "filename", "IC", "IC_universal", "length","consensus")  

datas=datas[,column_order]

#Add version
datas <- add_column(datas, version = "Version 1", .before = which(column_order=="type"))
datas=rename(datas, `type Version 1` = type)


# Version 2 motifs --------------------------------------------------------

#Some motifs are in Version 1 and Version 2, what is their type in Version 2
datas <- add_column(datas, `type Version 2` = NA, .before = which(names(datas)=="comment"))

datas1=datas

types=c("Batch_spacing_final_v2", "BatchOverlap_final_v2") #spacing composite
path="../../PWMs_final_version2/fromYimeng/"

#Dots were removed from the filenames, see ~/projects/TFBS/RProjects/TFBS/code/extract_motifs_from_excel_fromYimeng_version2.R

datas=tibble()

for(t in types){
  #t=types[1]
  data=as_tibble(do.call(rbind,strsplit(dir(paste0(path, t, "/")),split="_")))
  data$V8=gsub(".pfm", "", data$V8)
  data$symbol=paste(data$V1, data$V2,sep="_")
  data=data[,-c(1,2)]
  names(data)=c("ligand", "batch", "seed", "multinomial", "cycle", "short", "symbol")
  data=data %>% relocate(symbol, .before=ligand)
  data$study="fromYimeng"
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


append=TRUE #For writing the .scpd file

datas$ID=""
datas$filename_final=""
datas$IC=NA
datas$length=NA
datas$IC_universal=NA
datas$consensus=""
datas$`type Version 1`=NA
datas$`type Version 2`=NA

#Remove these from metadata, empty TFs
remove=c()

#Motifs that are already in version 1, remove these from datas
version1_motifs=c()

#Indexes of motifs in Version 2 that are in Version 1
match_indexes=c()

table(datas$type)
# Batch_spacing_final_v2  BatchOverlap_final_v2 
# 298                    971 

#All IDs are unique
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
    filename_final=paste0(gsub("_short.pfm", "", filename) ,".pfm")
    
    
    if(paste0(pwms,filename_final) %in% datas$filename_final){
      #There should not be such motifs
      print("Warning! Non-unique filenames and IDs")
      print(m)
    }
    
    ID=gsub(".pfm", "",filename_final)
    datas$ID[m]=ID
    
    if(datas$type[m]=="BatchOverlap_final_v2"){
      datas$`type Version 2`[m]="composite"
    }else if(datas$type[m] %in% c("Batch_spacing_final_v2")){
      datas$`type Version 2`[m]="spacing"
    }
    
    #Is the motif already in Version 1
    match_index=which(datas1$ID==ID)
    match_indexes=c(match_indexes, match_index)
    
    if(length(match_index)!=0){
      #Motif is already in Version 1
      version1_motifs=c(version1_motifs, m)
      #What is the motif type in Version 2
      if(datas$`type Version 2`[m]=="spacing"){
        datas1$`type Version 2`[match_index]="spacing"
      }else{
        #BatchOverlap_final_v2  
        datas1$`type Version 2`[match_index]="composite"
        
      }
      
      #Check that the pfms match
      version1_pfm=read.table(datas1$filename[match_index]) #.pfm is missing from the filename?
      version2_pfm=read.table(paste0("../../" ,datas$filename[m]))
      if( length(which(version1_pfm==version2_pfm))==nrow(version1_pfm)*ncol(version2_pfm) ){
        #Same motifs!
      }else{
        print(m)
        print("pfms of unique IDs are not the same!")
      }
      
      
      
    }else{
      #Motif is not in Version 1
      #Need to write the pfm files
      
      
      #write if to pwms, tab separated
      write.table(PWM,row.names = FALSE, col.names=FALSE, quote=FALSE,
                  file=paste0(pwms,filename_final) ,sep="\t")
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
                  file=paste0("../../PWMs_final_union/fromYimeng/Homo_sapiens_all", ".scpd"))
      append=TRUE
      
      write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
                  file=paste0("../../PWMs_final_union/fromYimeng/Homo_sapiens_all", ".scpd"))
      
      
      datas$filename_final[m]=paste0(pwms,filename_final)
    }
  }
}

#How many Version 1 motifs are also in Version 2?

table(datas1$`type Version 2`)
#missing          composite   spacing 
#149                    496    +    62 = 558

#Which indexes of datas1 are in Version 2
length(match_indexes) #558

#Which indexes of datas (Version2) are in datas1 (Version 1)
length(version1_motifs) #558

#Combine the two tables
#Do not include the motifs in Version 2 which are already in Version 1

datas=datas[-version1_motifs,]
datas$version="Version 2"

datas$filename=NULL
datas <- rename(datas, filename = filename_final)

#Add  family info

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
# "ZNF238" "C2H2 ZF"

commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB2")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB5")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB8")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="gntR")]="GntR-type HTH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="ZNF238")]="C2H2 ZF"

table(commondf_pairs2$DBD0, useNA="always") #There are no NAs
table(commondf_pairs2$DBD1, useNA="always") #There are no NAs


pair_families <- paste(commondf_pairs2$DBD0,commondf_pairs2$DBD1, sep = "_")

datas <- add_column(datas, Lambert2018_families = pair_families, .before = 4)

unique_Lambertfamilies <- unique(datas$Lambert2018_families)

datas$type=NULL

datas=datas[,names(datas1)]

#Concatenate data

datas_union=rbind(datas1, datas)



table(datas_union$version,useNA = "always")
# Version 1 Version 2      <NA> 
#  707       711         0 
table(datas_union$`type Version 1`,useNA = "always")
#both composite   spacing      <NA> 
#  2       501       204       711 (this many extra motifs in Version 2)
table(datas_union$`type Version 2`,useNA = "always")
# composite   spacing      <NA> 
#   971       298       149 (This many Version 1 motifs missing from Version 2)

#For motifs in both Version 1 and Version 2, for how many cases the composite/spacing info matches?
# The intersection of different versions is 558

table( !is.na(datas_union$`type Version 1`) & !is.na(datas_union$`type Version 2`) )
# FALSE  TRUE 
# 860    558 

#For how many motifs in the intersection the type information matches



table(!is.na(datas_union$`type Version 1`) & !is.na(datas_union$`type Version 2`) & datas_union$`type Version 2`==datas_union$`type Version 1`)
#FALSE  TRUE 
#878   540

# What are those for which the info does not match

type_info_index=which(!is.na(datas_union$`type Version 1`) & !is.na(datas_union$`type Version 2`))

table(datas_union$`type Version 2`[type_info_index]==datas_union$`type Version 1`[type_info_index])
#FALSE(this many motifs disagree in type between the versions)  TRUE 
#18   540

mismatch=type_info_index[which(datas_union$`type Version 2`[type_info_index]!=datas_union$`type Version 1`[type_info_index])]

tmp=datas_union[mismatch,]

#What about the computational definition of spacing and composite

capselex=readRDS(file="/Users/osmalama/projects/TFBS/RProjects/TFBS/RData/half_site_recognition_in_capselex_motifs.RDS")

fromYimeng=capselex %>% filter(study=="fromYimeng")

fromYimeng=fromYimeng[-grep(both_types[1], fromYimeng$ID)[1],]
fromYimeng=fromYimeng[-grep(both_types[2], fromYimeng$ID)[1],]

fromYimeng[grep(both_types[1], fromYimeng$ID)[1],"type"]="both"
fromYimeng[grep(both_types[2], fromYimeng$ID)[1],"type"]="both"

#These are all unique
fromYimeng$new_ID=paste0(fromYimeng$symbol, "_", fromYimeng$ligand, "_", fromYimeng$batch, "_", fromYimeng$seed, "_", fromYimeng$multinomial, "_", fromYimeng$cycle)

matches=match(datas_union$ID, fromYimeng$new_ID)
nona_ind=!is.na(matches)
head(datas_union$ID)
head(fromYimeng$new_ID[matches[nona_ind]])




datas_union <- add_column(datas_union, type_computation = NA, .before = which(names(datas_union)=="comment"))
datas_union$type_computation[nona_ind]=fromYimeng$type_computation[matches[nona_ind]]

#Which Version2 type and type_computation does not match?

nona_ind=which(!is.na(datas_union$`type Version 2`) & !is.na(datas_union$type_computation)  ) #558
  
table(datas_union$`type Version 2`[nona_ind]==datas_union$type_computation[nona_ind])
#FALSE (disagreement between Version 2 motifs and type_computation)  TRUE 
#129   429 

#Which Version1 type and type_computation does not match?

nona_ind=which(!is.na(datas_union$`type Version 1`) & !is.na(datas_union$type_computation)  ) #558
table(datas_union$`type Version 1`[nona_ind]==datas_union$type_computation[nona_ind])

#FALSE (disagreement between Version 2 motifs and type_computation)  TRUE 
#129   429 
write.table(datas_union, file="../../PWMs_final_union/fromYimeng/metadata.csv", row.names = FALSE,sep="\t")

saveRDS(datas_union, file="RData/fromYimeng_union.Rds")
datas_union=readRDS(file="RData/fromYimeng_union.Rds")

#Create data table for google drive

names(datas_union)

datas_union$clone=NULL
datas_union$family=NULL

datas_union <- rename(datas_union, family= Lambert2018_families)

datas_union$organism=NULL
datas_union$study=NULL
datas_union$experiment=NULL
datas_union$filename=NULL
datas_union$IC_universal=NULL
datas_union$short=NULL


datas_union$comment=""

datas_union$representative=""

datas_union <- rename(datas_union,  `type computation`=type_computation)

datas_union$`number of motif matches in GRCh38`=""
datas_union$`motif match threshold`=""

table(datas_union$`type Version 1`)
table(datas_union$`type Version 2`)

datas_union <- add_column(datas_union, `type computation overlap or gap` = "", .before = which(names(datas_union)=="comment"))

write.table(datas_union, file="../../PWMs_final_union/metadata_google_drive_supplement.tsv",quote=FALSE, sep="\t", row.names = FALSE)

stop()

#Missing from Ilya

library(readr)
absent <- read_csv("~/Downloads/abesnt.csv")
which(absent$ID %in% datas_union$ID)

metadata <- read_delim("~/projects/TFBS/PWMs_final/metadata.csv", 
                           delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

index=c()
for(i in 1:nrow(absent)){
    index[i]=grep(absent$ID[i], metadata$ID)
}

tmp=metadata[index,]


library("dplyr")
as_tibble(metadata) %>% filter(study=="fromYimeng") %>% count(type)



# 20230420             14
# pfm_composite_new   503
# pfm_spacing_new     192

# composite 503 and 206 spacing

# both composite   spacing 
# 2       501       204 

#composite 971   spacing 298
#971       298 


#datas_union$type_computation=NULL
#datas_union$clone=NULL
#datas_union$family=NULL


#write.table(datas_union, file="../../PWMs_final_union/metadata.tsv",quote=FALSE, sep="\t", row.names = FALSE)

nrow(datas_union) #1418


# Curated spacing and composite motifs ------------------------------------

library(readxl)
#Composites 391
curated_composites <- read_excel("~/projects/TFBS/curated_spacing_and_composite/TF_pairs_with_composite_binding_sites_filter.common20240117.xlsx")
curated_composites$curation="composite"
#Spacing 1613
curated_spacing <- read_excel("~/projects/TFBS/curated_spacing_and_composite/TF_pairs_showing_spacing_orientation_preference.filter_common_20240117.xlsx")
curated_spacing$curation="spacing"

curated_spacing=curated_spacing %>% mutate(ID2 = paste(HGNC, Barcode,Batch, Seed, paste0("m", Mul),paste0("c",Cyc),sep="_")) #These are unique

datas_union=datas_union %>%  mutate(ID2 = paste(symbol, ligand, batch, seed, multinomial,cycle,sep="_")) #Also unique
length(unique(datas_union$ID2)) #These are not unique 1418


notna_ind=which(!is.na(match( datas_union$ID2, curated_spacing$ID2))) #138
head(datas_union$ID2[notna_ind])
head(curated_spacing$ID2[match( datas_union$ID2, curated_spacing$ID2)[notna_ind]])

datas_union <- left_join(datas_union, curated_spacing, by = "ID2")
#New columns: "HGNC" "Length" "Barcode" "Batch" "Seed" "Mul" "Cyc" "DOMAIN" "1" "2" "curation.x (changed later)" 


curated_spacing=curated_spacing[-match( datas_union$ID2, curated_spacing$ID2)[notna_ind],] #left with 1475

#How many of the original spacing were curated as spacing
length(which(datas_union$`type Version 1`=="spacing" & datas_union$curation=="spacing")) # 127
length(which(datas_union$`type Version 2`=="spacing" & datas_union$curation=="spacing")) # 59
#How many of the original composites were curated as spacing
length(which(datas_union$`type Version 1`=="composites" & datas_union$curation=="spacing")) # 0
length(which(datas_union$`type Version 2`=="composites" & datas_union$curation=="spacing")) # 0
#For how many motifs, computational spacing and curated spacing agrees
length(which(datas_union$type_computation=="spacing" & datas_union$curation=="spacing")) # 28


#Try to match composites

#TF pairs that have ligand in their name
datas_union$ID_composite=paste0(datas_union$symbol,"_", datas_union$ligand) #These are not unique
length(unique(datas_union$ID_composite)) #1099

#Which are not unique
datas_union$ID[which(datas_union$ID_composite %in% names(which(table(datas_union$ID_composite)>1)))]

notna_ind=which(!is.na(match(curated_composites$TF_pairs, datas_union$ID_composite))) #204
head(curated_composites[notna_ind,])
head(datas_union$ID_composite[match(curated_composites$TF_pairs,datas_union$ID_composite)[notna_ind]])

#are these unique, NOT
length(unique(datas_union$ID_composite[match(curated_composites$TF_pairs,datas_union$ID_composite)[notna_ind]])) #200
#Which ones are not unique

not_unique_ind=which(datas_union$ID_composite[match(curated_composites$TF_pairs,datas_union$ID_composite)[notna_ind]] %in% names(which(table(datas_union$ID_composite[match(curated_composites$TF_pairs,datas_union$ID_composite)[notna_ind]])>1)))
datas_union$ID_composite[match(curated_composites$TF_pairs,datas_union$ID_composite)[notna_ind]][not_unique_ind]

#Remove the non-unique EN2_TBX20 is still remainin
curated_composites_cleaned=curated_composites[which(!(curated_composites$TF_pairs %in% unique(datas_union$ID_composite[match(curated_composites$TF_pairs,datas_union$ID_composite)[notna_ind]][not_unique_ind]))),]

notna_ind=which(!is.na(match(curated_composites_cleaned$TF_pairs, datas_union$ID_composite))) #196
head(curated_composites_cleaned[notna_ind,])
head(datas_union$ID_composite[match(curated_composites_cleaned$TF_pairs,datas_union$ID_composite)[notna_ind]])

length(unique(datas_union$ID_composite[match(curated_composites_cleaned$TF_pairs,datas_union$ID_composite)[notna_ind]])) #196


#Does anyone of these match to those already curated as spacing NO
#which( datas_union$ID_composite %in% datas_union$ID_composite[match(curated_composites$TF_pairs,datas_union$ID_composite)[notna_ind]] & datas_union$curation=="spacing") #0

curated_composites_cleaned$ID_composite=curated_composites_cleaned$TF_pairs

datas_union <- left_join(datas_union, curated_composites_cleaned, by = "ID_composite")
#New columns 
#"TF_pairs (TF_pairs.x)" "consensus.y" "p_value (p_value.x)" "comments (comments.x)" "Batch.y" "curation.y"

#There are some TF-TF-pairs which occur several times in datas_union but only once in curated_composites_cleaned, what are these and how many

matched_composites=paste0(datas_union$TF_pairs, datas_union$consensus.y, datas_union$p_value, datas_union$comments, datas_union$Batch.y, datas_union$curation.y)
matched_composites_names=names(which(table(matched_composites)>1))
matched_composites_names=matched_composites_names[-which(matched_composites_names=="NANANANANANA")]
datas_union$ID[which(matched_composites %in% matched_composites_names)] #There can be some motifs with both original composite and spacing type


#Remove those that were already matched from curated_composites
#tmp=curated_composites_cleaned #383
#322
curated_composites_cleaned=curated_composites_cleaned[-which(paste0(curated_composites_cleaned$TF_pairs, curated_composites_cleaned$consensus, 
                                                                    curated_composites_cleaned$p_value, curated_composites_cleaned$comments, curated_composites_cleaned$Batch, curated_composites_cleaned$curation) %in% matched_composites_names),]

#Match the symbols





#TF pairs have ligand in their name
notna_ind=which(!is.na(match(curated_composites_cleaned$TF_pairs, datas_union$symbol))) #162
head(curated_composites_cleaned[notna_ind,])
head(datas_union$symbol[match(curated_composites_cleaned$TF_pairs,datas_union$symbol)[notna_ind]])

length(unique(curated_composites_cleaned$TF_pairs[notna_ind])) #160

names(which(table(curated_composites_cleaned$TF_pairs[notna_ind])>1))
names(which(table(datas_union$symbol[match(curated_composites_cleaned$TF_pairs,datas_union$symbol)[notna_ind]])>1))

#"ELK1_TBX21"  "PAX2_NKX6-1"

#Remove the non-unique

not_unique_ind=which(datas_union$symbol[match(curated_composites_cleaned$TF_pairs,datas_union$symbol)[notna_ind]] %in% names(which(table(datas_union$symbol[match(curated_composites_cleaned$TF_pairs,datas_union$symbol)[notna_ind]])>1)))
curated_composites_cleaned2=curated_composites_cleaned[which(!(curated_composites_cleaned$TF_pairs %in% unique(datas_union$symbol[match(curated_composites_cleaned$TF_pairs,datas_union$symbol)[notna_ind]][not_unique_ind]))),]

curated_composites_cleaned2$symbol=curated_composites_cleaned2$TF_pairs

#Are there some TF_pairs for which we have already obtained all info, remove these from curated_composites_cleaned2

for(pair in curated_composites_cleaned2$symbol[is.na(curated_composites_cleaned2$consensus)]){
  #pair=curated_composites_cleaned2$symbol[is.na(curated_composites_cleaned2$consensus)][1]
  if(length(which(datas_union$symbol %in% pair))!=0){
    #print(pair)   
    if(names(table(!is.na(datas_union$TF_pairs[ which(datas_union$symbol==pair)]))) %in% c("TRUE")
       
       #!is.na(datas_union$TF_pairs[ which(datas_union$symbol==pair)])
    ){
      print(pair)
    }
  }
  
}

datas_union %>% filter(symbol=="HOXA6_TBX4") %>% select(TF_pairs, consensus.y, p_value, comments, Batch.y, curation.y)
datas_union %>% filter(symbol=="HOXB7_TBX4")  %>% select(TF_pairs, consensus.y, p_value, comments, Batch.y, curation.y)
datas_union %>% filter(symbol=="VSX1_TBX4")  %>% select(TF_pairs, consensus.y, p_value, comments, Batch.y, curation.y)
datas_union %>% filter(symbol=="HOXA4_TBX20")  %>% select(TF_pairs, consensus.y, p_value, comments, Batch.y, curation.y)
datas_union %>% filter(symbol=="HOXB5_TBX20")  %>% select(TF_pairs, consensus.y, p_value, comments, Batch.y, curation.y)
datas_union %>% filter(symbol=="HOXD10_TBX20")  %>% select(TF_pairs, consensus.y, p_value, comments, Batch.y, curation.y)

curated_composites_cleaned2 =curated_composites_cleaned2[-which(curated_composites_cleaned2$TF_pairs %in% c("HOXA6_TBX4", "HOXB7_TBX4", "VSX1_TBX4", "HOXA4_TBX20", "HOXB5_TBX20", "HOXD10_TBX20", "PAX5_EVX2","TBX4_EVX2")),]



datas_union <- left_join(datas_union, curated_composites_cleaned2, by = "symbol")




#New columns: "TF_pairs.y"  "consensus" "p_value.y" "comments.y"  "Batch" "curation"  "ID_composite.y"      


#Is there some motifs that are not found in the curated set


which(!(is.na(datas_union$TF_pairs.y)))

datas_union[which(!(is.na(datas_union$TF_pairs.y))) , c("TF_pairs.x", "consensus.y", "p_value.x",  "comments.x", "Batch.y", "curation.y")] = datas_union[which(!(is.na(datas_union$TF_pairs.y))), c("TF_pairs.y", "consensus", "p_value.y",
                                                                                                                                                                                                "comments.y","Batch", "curation")   ]

datas_union=datas_union[, -which(names(datas_union) %in% c("TF_pairs.y", "consensus", "p_value.y", "comments.y", "Batch", "curation","ID_composite.y" ))]

datas_union <- datas_union %>%
  rename(
    consensus = consensus.x,
    HGNC_spacing=HGNC,
    Length_spacing=Length,
    Barcode_spacing=Barcode,
    Batch_spacing=Batch.x,                    
    Seed_spacing=Seed,
    Mul_spacing=Mul,
    Cyc_spacing=Cyc,
    DOMAIN_spacing=DOMAIN,
    one_spacing="1",
    two_spacing="2",                           
    curation_spacing=curation.x,
    ID_composite=ID_composite.x,
    TF_pairs_composite=TF_pairs.x,
    consensus_composite=consensus.y,
    p_values_composite=p_value.x,
    comments_composite=comments.x,
    Batch_composite=Batch.y,
    curation_composite=curation.y                 
  )






test=datas_union %>% dplyr::select(ID, symbol, Lambert2018_families, ligand, batch, seed, multinomial, filename, consensus,version,`type Version 1`,       `type Version 2`,       type_computation,  curation_spacing, curation_composite, 
                            ID2,HGNC_spacing, Length_spacing, Barcode_spacing, Batch_spacing,              
                           Seed_spacing,                Mul_spacing,                 Cyc_spacing,               DOMAIN_spacing,              one_spacing,                 two_spacing,                
                                      ID_composite,                TF_pairs_composite,        consensus_composite,         p_values_composite,          comments_composite,         
                           Batch_composite    )

library(openxlsx)

# Write the tibble to an Excel file
write.xlsx(test, "../../PWMs_final_union/final_curation.xlsx")
########### Old stuff below #########################


