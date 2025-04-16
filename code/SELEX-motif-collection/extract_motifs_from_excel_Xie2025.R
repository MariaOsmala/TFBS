

library(tidyverse)
library(readxl)
library(universalmotif)
library(ggplot2)
rm(list=ls())
source("read_excel_tables.R")
source("paths_to_folders.R", echo=FALSE)
source("compute_kl_divergence.R", echo=FALSE)
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


all_table <- read_excel("../../Data/SELEX-motif-collection/41586_2025_8844_MOESM5_ESM.xlsx", 
                        sheet="Supplementary Table S3", col_names=FALSE, skip=13, col_types=c(rep("text",9), rep("numeric",20)))






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
PWMs_metadata$kld_between_revcomp=NA #Kullback-Leibler divergence between the motif and its reverse complement. Gives idea whether the motif is palindromic.


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
                                                                              "representative", "short", "type", "filename","ID", "IC","length", "consensus", "kld_between_revcomp"))], collapse="_"),
                    ".pfm")
    ID=paste0(PWMs_metadata[m,
                            -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "experiment","representative", "short", "type", "filename",
                                                                 "ID", "IC", "length", "consensus", "kld_between_revcomp" ))], collapse="_")
    
    if(filename %in% PWMs_metadata$filename){
      print("Warning! Non-unique filenames and IDs")
      print(ID)
      filename=paste0(pfms_tab_path,"/", 
                      paste0(PWMs_metadata[m,
                                           -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "experiment",
                                                                                "representative", "short", "type", "filename","ID", "IC", "length", "consensus", "kld_between_revcomp"))], collapse="_"),
                      "_v2",".pfm")
      ID=paste0( paste0(PWMs_metadata[m,
                                      -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "experiment","representative", "short", "type", "filename",
                                                                           "ID", "IC", "length", "consensus","kld_between_revcomp" ))], collapse="_"), "_v2")
      
      
      
    }
    
    PWMs_metadata$filename[m]=filename
    PWMs_metadata$ID[m]=ID
    
    
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,file=filename, sep="\t")
    
    pcm=as.matrix( read.csv(filename, header=FALSE, sep="\t") )
    colnames(pcm)=NULL
    rownames(pcm)=c("A", "C", "G", "T")
    
    # Make reverse complement
    # Reverse the matrix = flip the matrix horizontally
    pcm_rev=pcm[,ncol(pcm):1]
    
    #Substitute nucleotides with complements = flip the orders of A&T and C&G
    
    pcm_rev_comp=pcm_rev
    pcm_rev_comp["A", ]=pcm_rev["T",]
    pcm_rev_comp["T", ]=pcm_rev["A",]
    pcm_rev_comp["C", ]=pcm_rev["G",]
    pcm_rev_comp["G", ]=pcm_rev["C",]
    
    
    
    pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),name=PWMs_metadata$symbol[m],
                                profileMatrix=pcm
    )
    
    pfm_rev_comp <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),name=PWMs_metadata$symbol[m],
                                         profileMatrix=pcm_rev_comp
    )
    
    
    pwm_class <- TFBSTools::toPWM(pfm, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
    
    pwm_class_revcomp <- TFBSTools::toPWM(pfm_rev_comp, type="prob", 
                                          pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
    
    #Compute the KL-divergence between the original and the reverse complement
    
    # Compute KL divergence, it does not matter in which order you give the matrices to the function
    kl_divergence <- as.numeric(compute_kl_divergence(pwm_class@profileMatrix, pwm_class_revcomp@profileMatrix))
    
    PWMs_metadata[m, "kld_between_revcomp"]=kl_divergence
    
    
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
    
    #Write reverse complements
    
    write.table(pcm_rev_comp, file=paste0(rev_comp_pfms_space_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
    write.table(pcm_rev_comp, file=paste0(rev_comp_pfms_tab_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
    
    
    
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
    
    
    pdf(paste0(logos_pdf_ic,"/",ID,".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto" ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    dev.off()
    
    
    #Reverse complements 
    
    png(paste0(revcomp_logos_png_prob,"/",ID,".png"), res=600,  width = 2500/4*ncol(pwm_class_revcomp@profileMatrix), height = 2500)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    pdf(paste0(revcomp_logos_pdf_prob,"/",ID,".pdf"), width=ncol(pwm_class_revcomp@profileMatrix), height = 4)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto" ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    png(paste0(revcomp_logos_png_ic,"/",ID,".png"), res=600,  width = 2500/4*ncol(pwm_class_revcomp@profileMatrix), height = 2500)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    pdf(paste0(revcomp_logos_pdf_ic,"/",ID,".pdf"),  width=ncol(pwm_class_revcomp@profileMatrix), height = 4)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
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
                "type", "comment",  "filename", "IC",  "length","consensus", "kld_between_revcomp")  

PWMs_metadata=PWMs_metadata[,column_order]



# Novel HT-SELEX monomers----------------------------------------------------------

PWMs=split_df(all_table)[[3]]

monomer_metadata=do.call(rbind,lapply(
  seq(1,nrow(PWMs),5), 
  function(i) PWMs[i,-1]
)
)


colnames(monomer_metadata)=c("symbol",	"ligand",	"batch", "seed",	"multinomial",
                          "cycle","type","representative")


monomer_metadata=monomer_metadata[, c("symbol",	"ligand",	"batch", "seed",	"multinomial",
                                "cycle","type","representative")]


monomer_metadata$clone=NA
monomer_metadata$study="Xie2025"
monomer_metadata$experiment="HT-SELEX"
monomer_metadata$short=NA
monomer_metadata$family=NA
monomer_metadata$filename=NA
monomer_metadata$comment=NA
monomer_metadata$organism="Homo_sapiens"



monomer_metadata=select(monomer_metadata,c("symbol",	"clone","family", "organism",	"study","experiment",
                                     "ligand",	"batch", "seed",	"multinomial",
                                     "cycle","representative", "short", "type","comment", "filename"))


monomer_metadata$ID=""
monomer_metadata$IC=NA
monomer_metadata$length=NA
monomer_metadata$consensus=""
monomer_metadata$kld_between_revcomp=NA #Kullback-Leibler divergence between the motif and its reverse complement. Gives idea whether the motif is palindromic.


PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )



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
                    paste0(monomer_metadata[m,
                                         -which(colnames(monomer_metadata)%in% c("clone", "family","organism", "study","comment", "experiment",
                                                                              "representative", "short", "type", "filename","ID", "IC","length", "consensus", "kld_between_revcomp"))], collapse="_"),
                    ".pfm")
    ID=paste0(monomer_metadata[m,
                            -which(colnames(monomer_metadata)%in% c("clone", "family","organism", "study","comment", "experiment","representative", "short", "type", "filename",
                                                                 "ID", "IC", "length", "consensus", "kld_between_revcomp" ))], collapse="_")
    
    if(filename %in% monomer_metadata$filename){
      print("Warning! Non-unique filenames and IDs")
      print(ID)
      filename=paste0(pfms_tab_path,"/", 
                      paste0(monomer_metadata[m,
                                           -which(colnames(monomer_metadata)%in% c("clone", "family","organism", "study","comment", "experiment",
                                                                                "representative", "short", "type", "filename","ID", "IC", "length", "consensus", "kld_between_revcomp"))], collapse="_"),
                      "_v2",".pfm")
      ID=paste0( paste0(monomer_metadata[m,
                                      -which(colnames(monomer_metadata)%in% c("clone", "family","organism", "study","comment", "experiment","representative", "short", "type", "filename",
                                                                           "ID", "IC", "length", "consensus","kld_between_revcomp" ))], collapse="_"), "_v2")
      
      
      
    }
    
    monomer_metadata$filename[m]=filename
    monomer_metadata$ID[m]=ID
    
    
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,file=filename, sep="\t")
    
    pcm=as.matrix( read.csv(filename, header=FALSE, sep="\t") )
    colnames(pcm)=NULL
    rownames(pcm)=c("A", "C", "G", "T")
    
    # Make reverse complement
    # Reverse the matrix = flip the matrix horizontally
    pcm_rev=pcm[,ncol(pcm):1]
    
    #Substitute nucleotides with complements = flip the orders of A&T and C&G
    
    pcm_rev_comp=pcm_rev
    pcm_rev_comp["A", ]=pcm_rev["T",]
    pcm_rev_comp["T", ]=pcm_rev["A",]
    pcm_rev_comp["C", ]=pcm_rev["G",]
    pcm_rev_comp["G", ]=pcm_rev["C",]
    
    
    
    pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),name=monomer_metadata$symbol[m],
                                profileMatrix=pcm
    )
    
    pfm_rev_comp <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),name=monomer_metadata$symbol[m],
                                         profileMatrix=pcm_rev_comp
    )
    
    
    pwm_class <- TFBSTools::toPWM(pfm, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
    
    pwm_class_revcomp <- TFBSTools::toPWM(pfm_rev_comp, type="prob", 
                                          pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
    
    #Compute the KL-divergence between the original and the reverse complement
    
    # Compute KL divergence, it does not matter in which order you give the matrices to the function
    kl_divergence <- as.numeric(compute_kl_divergence(pwm_class@profileMatrix, pwm_class_revcomp@profileMatrix))
    
    monomer_metadata[m, "kld_between_revcomp"]=kl_divergence
    
    
    icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
    monomer_metadata[m, "IC"]=sum(rowSums(icm)) 
    
    #motif length
    monomer_metadata[m, "length"]=length(pfm)
    
    motif=universalmotif::read_matrix(file=filename, sep="\t", header=FALSE)
    motif@name=ID
    
    monomer_metadata[m, "consensus"]=motif@consensus
    
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
    
    #Write reverse complements
    
    write.table(pcm_rev_comp, file=paste0(rev_comp_pfms_space_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
    write.table(pcm_rev_comp, file=paste0(rev_comp_pfms_tab_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
    
    
    
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
    
    
    pdf(paste0(logos_pdf_ic,"/",ID,".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto" ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    dev.off()
    
    
    #Reverse complements 
    
    png(paste0(revcomp_logos_png_prob,"/",ID,".png"), res=600,  width = 2500/4*ncol(pwm_class_revcomp@profileMatrix), height = 2500)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    pdf(paste0(revcomp_logos_pdf_prob,"/",ID,".pdf"), width=ncol(pwm_class_revcomp@profileMatrix), height = 4)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto" ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    png(paste0(revcomp_logos_png_ic,"/",ID,".png"), res=600,  width = 2500/4*ncol(pwm_class_revcomp@profileMatrix), height = 2500)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    pdf(paste0(revcomp_logos_pdf_ic,"/",ID,".pdf"),  width=ncol(pwm_class_revcomp@profileMatrix), height = 4)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    
  }
}


 
# Merge two data frames
monomer_metadata <- monomer_metadata %>%
  left_join(select(TF_family_mappings, "symbol", "Gene Information/ID","DBD" ), by = 'symbol')

monomer_metadata <- add_column(monomer_metadata, Lambert2018_families = monomer_metadata$DBD, .before = 4)


monomer_metadata <- monomer_metadata %>% #Save this info only for human monomers
  dplyr::rename("Human_Ensemble_ID"= "Gene Information/ID") %>%
  relocate("Human_Ensemble_ID", .after = "symbol")

monomer_metadata$DBD=NULL







monomer_metadata[which(is.na(monomer_metadata$Lambert2018_families)),]
#THHEX, do not know the family of this 



PWMs_metadata$Human_Ensemble_ID=NA
PWMs_metadata<- PWMs_metadata %>%
  relocate("Human_Ensemble_ID", .after = "symbol")

names(PWMs_metadata)[which(!is.element(names(PWMs_metadata),names(monomer_metadata)))]

monomer_metadata=monomer_metadata[,names(PWMs_metadata)]

PWMs_metadata=rbind(PWMs_metadata, monomer_metadata)

write.table(PWMs_metadata, file="../../Data/SELEX-motif-collection/Xie2025_metadata.csv", row.names = FALSE,sep="\t")

