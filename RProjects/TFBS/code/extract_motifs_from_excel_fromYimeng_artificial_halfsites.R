
library(tidyverse)
library(universalmotif)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)

#I checked that the following TFs have a spacing motif but the monomer looks homodimer. 
#It would be useful to have the true monomeric half-site motifs for these. 
#For SOX11 and PROX1, I don't know if the half-sites in the homodimer are equal.

# IRF8
# MAFF
# NFATC1
# FOXD2
# SOX11 ?
# TGIF2LX
# ONECUT1
# PROX1 ?
# TFAP2C

#All the TFs in your list prefer to bind to DNA as homodimer in HT-SELEX 
#and don't properly bind to DNA as a monomer. 

#Anyhow I generated half-site motif for these TFs or their paralogs 
#(FOXD3 PWM for FOXD2, IRF6 PWM for IRF8) 
#except PROX1 and TFAP2C, for which the half-site motifs are impossible to generate.  
#I am not sure the PROX1 motif we have in our hand is classified as monomer or homodimer. 

#datas$symbol "FOXD3"   "IRF6"    "MAFF"    "NFATC1"  "ONECUT1" "SOX11"   "TGIF2LX"

rm(list=ls())

#Version2.2

#HT-SELEX
length(dir("/Users/osmalama/projects/TFBS/PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/")) #7

path="../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/"

datas=tibble()


data=as_tibble(do.call(rbind,strsplit(dir(path)[grep(".pfm",dir(path))],split="_")))
data$V7=gsub(".pfm", "", data$V7)
names(data)=c("symbol","ligand", "batch", "seed", "multinomial", "cycle", "short")
data$study="artificial_halfsites"
data$type=NA
data$organism="Homo_sapiens"
data$clone=NA
data$family=NA
data$experiment="HT-SELEX"
data$representative=NA
data$comment=NA

data$filename=paste0("PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/", dir(path)[grep(".pfm",dir(path))])
data$ID=gsub(".pfm", "",do.call(rbind, strsplit(data$filename, "/"))[,4])
data=data %>% relocate(ID, .before=symbol)
datas=rbind(datas, data)

datas$symbol


datas=datas[, c("ID","symbol",	"clone","family",	"organism","study","experiment",
                "ligand",	"batch", "seed",	"multinomial",
                "cycle","short","representative", "type","comment","filename")]


pwms=paste0("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/pwms/")
pwms_space=paste0("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/pwms_space/")
transfac_path=paste0("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/pwms/transfac/")
dir.create(pwms, recursive=TRUE)
dir.create(pwms_space, recursive=TRUE)
dir.create(transfac_path, recursive=TRUE)

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
      filename=paste0(strsplit(f, "/")[[1]][-c(1,2,3)],collapse="/")  
    }else{
    print("error")
    }
    
    #add extension
    if(length(grep("_short.pfm", filename))!=0){
      filename_final=paste0(gsub("_short.pfm", "", filename), ".pfm") 
    }else{
      filename_final=filename
    }
    
    
  
    ID=gsub(".pfm", "",filename_final)
    datas$ID[m]=ID
        
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
                file=paste0("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/Homo_sapiens_all", ".scpd"))
    append=TRUE
    
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
                file=paste0("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/Homo_sapiens_all", ".scpd"))
    
    datas$filename_final[m]=paste0(pwms,filename_final)
  }
}



#Add family mappings

TF_family_mappings <- read_delim("~/projects/motif-clustering-Viestra-private/HumanTFs-ccbr-TableS1.csv", delim=",")


pairs <- as.data.frame(do.call(rbind, strsplit(datas$symbol, split="_")))
colnames(pairs)=c("symbol")
pairs$order=seq(nrow(pairs))

commondf_pairs <- merge(pairs, TF_family_mappings, by = 'symbol', all.x = TRUE)[,c(1,2,3,4,5)]
commondf_pairs <- commondf_pairs[order(commondf_pairs$order), ]



datas <- add_column(datas, Lambert2018_families = commondf_pairs$DBD, .before = 4)
datas <- rename(datas, filename = filename_final)

column_order=c( "ID", "symbol", "clone","family", "Lambert2018_families", "organism", "study","experiment",          
                "ligand", "batch","seed", "multinomial","cycle","representative",  "short",               
                "type", "comment",  "filename", "IC", "IC_universal", "length","consensus")  

datas=datas[,column_order]



write.table(datas, file="../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/metadata.csv", row.names = FALSE,sep="\t")

