
#library(dplyr)
#library(tidyr)
#library(purrr)
#library(stringr)
#library(readr)
library(tidyverse)
library(readxl)
library(universalmotif)
library(ggplot2)
rm(list=ls())
source("read_excel_tables.R")
source("paths_to_folders.R", echo=FALSE)

# Table S2 PWM models																									
# This table contains background subtracted PWMs for all of the TF pairs and individual TFs analyzed. 
# First line for each factor contains information as follows:																									

# Base;	A, C, G or T																								
# symbol(s)	HGNC symbol(s) (gene names); Either for reference models that were generated for individual TFs in this study or for the co-operative pairs, in which case the TF names are separated with an underscore "_". First of the TFs is TF1 that was expressed as a His tagged clone and was used in first of the affinity separations and the second one is the SBP tagged TF2																								
# ligand sequence	Name of primer used to generate selection ligand; number before N indicates length of random sequence																								
# batch;	batch of HT-SELEX or CAP-SELEX experiment																								
# seed;	Sequence used to identify subsequences to be included in the model																								
# multinomial;	Hamming distance to seed used to identify subsequences to be included in the model																								
# cycle	SELEX, cycle used to construct the model; previous cycle was used as background unless indicated otherwise by a letter "b" followed by the used background cycle, "u" letter indicates that only the first instance of identical sequencing reads were used in the construction of the PWM																								
# Interaction type;	composite or spacing 																								
# Matrix is one of the representative PWMs; value yes or no																								


all_table <- read_excel("../../Data/SELEX-motif-collection/Xie et al. 2025 Supplementary_Tables.xlsx", 
                        sheet="Table S2PWMs", col_names=FALSE, skip=13, col_types=c(rep("text",10), rep("numeric",19)))






PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(
  seq(1,nrow(PWMs),5), 
  function(i) PWMs[i,-1]
  )
)


colnames(PWMs_metadata)=c("symbol",	"ligand",	"batch", "seed",	"multinomial",
                          "cycle","type","representative")


PWMs_metadata=PWMs_metadata[, c("symbol",	"ligand",	"batch", "seed",	"multinomial",
                                "cycle","type","representative")]


PWMs_metadata$clone=NA
PWMs_metadata$study="Xie2025"
PWMs_metadata$experiment="CAP-SELEX"
PWMs_metadata$short=NA
PWMs_metadata$family=NA
PWMs_metadata$filename=NA
PWMs_metadata$comment=NA
PWMs_metadata$organism="Homo_sapiens"

PWMs_metadata$representative=toupper(PWMs_metadata$representative)

PWMs_metadata=select(PWMs_metadata,c("symbol",	"clone","family", "organism",	"study","experiment",
                                     "ligand",	"batch", "seed",	"multinomial",
                                     "cycle","representative", "short", "type","comment", "filename"))


PWMs_metadata$ID=""
PWMs_metadata$IC=NA
PWMs_metadata$length=NA
PWMs_metadata$consensus=""


PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )


append=FALSE


remove=c()

for(m in 1:length(PWMs_list)){
  
  numeric_matrix <- matrix(as.numeric(as.vector(as.matrix(PWMs_list[[m ]][,-1]))), nrow = nrow(as.matrix(PWMs_list[[m ]][,-1])))
  
  if( ncol(PWMs_list[[m ]]) ==1 ) {
    remove=c(remove, m)
    print(m)
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWMs_list[[ m ]][,-1]))) ==0 ){
    remove=c(remove, m)    
    print(m)
    
  }else{
    #There are zero columns in the matrix, remove them
    if(length(which(colSums(numeric_matrix)==0))>0){
      # Assume 'char_matrix' is your character matrix
      numeric_matrix <- matrix(as.numeric(as.vector(as.matrix(PWMs_list[[m ]][,-1]))), nrow = nrow(as.matrix(PWMs_list[[m ]][,-1])))
      print(m)
      #print(PWMs_list[[m ]])
      #Remove the zero column, it is always the last
      PWMs_list[[m ]]=PWMs_list[[m ]][,-ncol(PWMs_list[[m ]])]
      
    }
    
    
    
    
    # Write .pfm tab-separated
    filename=paste0(pfms_tab_path,"/", 
                    paste0(PWMs_metadata[m,
                                         -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "experiment",
                                                                              "representative", "short", "type", "filename","ID", "IC","length", "consensus"))], collapse="_"),
                    ".pfm")
    ID=paste0(PWMs_metadata[m,
                            -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "experiment","representative", "short", "type", "filename",
                                                                 "ID", "IC", "length", "consensus" ))], collapse="_")
    
    if(filename %in% PWMs_metadata$filename){
      print("Warning! Non-unique filenames and IDs")
      print(ID)
      filename=paste0(pfms_tab_path,"/", 
                      paste0(PWMs_metadata[m,
                                           -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "experiment",
                                                                                "representative", "short", "type", "filename","ID", "IC", "length", "consensus"))], collapse="_"),
                      "_v2",".pfm")
      ID=paste0( paste0(PWMs_metadata[m,
                                      -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "experiment","representative", "short", "type", "filename",
                                                                           "ID", "IC", "length", "consensus" ))], collapse="_"), "_v2")
      
      
      
    }
    
    PWMs_metadata$filename[m]=filename
    PWMs_metadata$ID[m]=ID
    
    
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,file=filename, sep="\t")
    
    pcm=as.matrix( read.csv(filename, header=FALSE, sep="\t") )
    colnames(pcm)=NULL
    rownames(pcm)=c("A", "C", "G", "T")
    
    
    
    pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),name=PWMs_metadata$symbol[m],
                                profileMatrix=pcm
    )
    
    
    pwm_class <- TFBSTools::toPWM(pfm, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
    
    
    icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
    PWMs_metadata[m, "IC"]=sum(rowSums(icm)) 
    
    #motif length
    PWMs_metadata[m, "length"]=length(pfm)
    
    motif=universalmotif::read_matrix(file=filename, sep="\t", header=FALSE)
    motif@name=ID
    
    PWMs_metadata[m, "consensus"]=motif@consensus
    
    transfac=paste0(pfms_transfac_path,"/", ID,".pfm")
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE, 
                file=paste0(pfms_space_path,"/", ID, ".pfm"), sep=" ")
    
    PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
                append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(pfms_scpd,"/","/Xie2025.scpd"))
    append=TRUE
    
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
                file=paste0(pfms_scpd,"/","/Xie2025.scpd"))
    
    
    #Write also pwms 
    
    write.table(pwm_class@profileMatrix, file=paste0(pwms_space_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
    write.table(pwm_class@profileMatrix, file=paste0(pwms_tab_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
    
    
    # Draw logos --------------------------------------------------------------
    
    
    
    png(paste0(logos_png_prob,"/",ID,".png"), res=600,  width = 2500/4*ncol(pwm_class@profileMatrix), height = 2500)
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    pdf(paste0(logos_pdf_prob,"/",ID,".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto" ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    dev.off()
    
    #Information content matrices
    
    
    png(paste0(logos_png_ic,"/",ID,".png"), res=600,  width = 2500/4*ncol(pwm_class@profileMatrix), height = 2500)
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    pdf(paste0(logos_pdf_ic,"/","ID",".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto" ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    dev.off()
    
    
  }
}

#Add family mappings

TF_family_mappings <- read_delim("../../Data/SELEX-motif-collection/HumanTFs-ccbr-TableS1.csv", delim=",")


pairs <- as.data.frame(do.call(rbind, strsplit(PWMs_metadata$symbol, split="_")))
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

#"THHEX" This is monomer?
#"BHLHB2" bHLH
#"BHLHB8" bHLH
#"FOXO3A" Forkhead
# "gntR"  GntR-type HTH domain 
# "ZNF238" "C2H2 ZF"
#"BHLHB5" bHLH

#These were obtained from UniProt
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB2")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB5")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="BHLHB8")]="bHLH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="FOXO3A")]="Forkhead"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="gntR")]="GntR-type HTH"
commondf_pairs2$DBD0[which(commondf_pairs2$symbol0=="ZNF238")]="C2H2 ZF"


table(commondf_pairs2$DBD0, useNA="always") #There is no NA
table(commondf_pairs2$DBD1, useNA="always") #There is no NA


pair_families <- paste(commondf_pairs2$DBD0,commondf_pairs2$DBD1, sep = "_")


PWMs_metadata <- add_column(PWMs_metadata, Lambert2018_families = pair_families, .before = 4)




column_order=c( "ID", "symbol", "clone","family", "Lambert2018_families", "organism", "study","experiment",          
                "ligand", "batch","seed", "multinomial","cycle","representative",  "short",               
                "type", "comment",  "filename", "IC",  "length","consensus")  

PWMs_metadata=PWMs_metadata[,column_order]



# Novel HT-SELEX ----------------------------------------------------------

path="../../Data/SELEX-motif-collection/Xie2025-HT-SELEX/" #12 motifs

datas=tibble()
data=as_tibble(do.call(rbind,strsplit(dir(path),split="_")))
data$V7=gsub(".pfm", "", data$V7)
names(data)=c("symbol","ligand", "batch", "seed", "multinomial", "cycle", "short")
data$study="Xie2025"
data$type=NA
data$organism="Homo_sapiens"
data$clone=NA
data$family=NA
data$experiment="HT-SELEX"
data$representative=NA
data$comment=NA
data$filename=paste0(path, dir(paste0(path)))
data$ID=gsub(".pfm", "",do.call(rbind, strsplit(data$filename, "/"))[,6])
data=data %>% relocate(ID, .before=symbol)
datas=rbind(datas, data)




datas=datas[, c("ID","symbol",	"clone","family",	"organism","study","experiment",
                "ligand",	"batch", "seed",	"multinomial",
                "cycle","short","representative", "type","comment","filename")]


datas <- datas %>%
  left_join(select(TF_family_mappings, "symbol", "Gene Information/ID","DBD" ), by = 'symbol')

datas <- add_column(datas, Lambert2018_families = datas$DBD, .before = 4)

#Rename Gene Information/ID column and relocate column 


datas<- datas %>% #Save this info only for human monomers
  dplyr::rename("Human_Ensemble_ID"="Gene Information/ID") %>%
  relocate("Human_Ensemble_ID", .after = "symbol")

datas$DBD=NULL


datas[which(is.na(datas$Lambert2018_families)),]
#THHEX, do not know the family of this 



PWMs_metadata$Human_Ensemble_ID=NA
PWMs_metadata<- PWMs_metadata %>%
  relocate("Human_Ensemble_ID", .after = "symbol")

names(PWMs_metadata)[which(!is.element(names(PWMs_metadata),names(datas)))]

append=TRUE

datas$IC=NA
datas$length=NA
datas$consensus=""
#Remove these from metadata, empty TFs
remove=c()


for(m in 1:nrow(datas)){
  print(m)
  f=datas$filename[m]
  PWM=as.matrix(read_tsv(f, col_names=FALSE))
  #Empty matrix
  if( ncol(PWM) ==1 ) {
    remove=c(remove, m)
    print(paste0("Empty matrix: ", f))
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWM[,-1]))) ==0 ){
    remove=c(remove, m)    
    print(paste0("Zero matrix: ",f))
  }else{
  
    
    filename=paste0(strsplit(f, "/")[[1]][-c(1,2,3,4,5)],collapse="/")  
    
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
                file=paste0(pfms_tab_path,"/",filename_final) ,sep="\t")
    
    
  
    write.table(PWM,row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(pfms_space_path,"/",filename_final) ,sep=" ")
    
    
    pcm=as.matrix( read.csv(paste0(pfms_tab_path,"/",filename_final), header=FALSE, sep="\t") )
    colnames(pcm)=NULL
    rownames(pcm)=c("A", "C", "G", "T")
    
    pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25), name=datas$symbol[m],
                                profileMatrix=pcm
    )
    
    icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
   
    datas$IC[m]=sum(rowSums(icm)) #6.77892 # universalmotif computes this wrong?
    
    #motif length
    datas$length[m]=length(pfm)
    
    motif=universalmotif::read_matrix(file=paste0(pfms_tab_path,"/",filename_final), sep="\t", header=FALSE)
    
   
    datas[m, "consensus"]=motif@consensus
    
    
    #Write transfac format
    
    transfac=paste0(pfms_transfac_path,"/",filename_final)
    
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
    
    #Write .scpd format, all into a single file
    PWM=as.matrix(PWM, dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
                append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(pfms_scpd,"/","/Xie2025.scpd"))
    append=TRUE
    
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
                file=paste0(pfms_scpd,"/","/Xie2025.scpd"))
    
    
    datas$filename[m]=paste0(pfms_tab_path,"/",filename_final)
    
  
  }
}







column_order=c( "ID", "symbol", "Human_Ensemble_ID","clone","family", "Lambert2018_families", "organism", "study","experiment",          
                "ligand", "batch","seed", "multinomial","cycle","representative",  "short",               
                "type", "comment",  "filename", "IC", "length","consensus")  
datas=datas[,column_order]

PWMs_metadata=rbind(PWMs_metadata, datas)



write.table(PWMs_metadata, file="../../Data/SELEX-motif-collection/Xie2025_metadata.csv", row.names = FALSE,sep="\t")

