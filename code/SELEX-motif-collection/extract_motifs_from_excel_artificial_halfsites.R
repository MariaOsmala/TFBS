
library(tidyverse)
library(readxl)
library(universalmotif)
library(ggplot2)
rm(list=ls())

source("compute_kl_divergence.R", echo=FALSE)


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

#All the TFs in this list prefer to bind to DNA as homodimer in HT-SELEX  and don't properly bind to DNA as a monomer. 

#Anyways, we generated half-site motif for these TFs or their paralogs (FOXD3 PWM for FOXD2, IRF6 PWM for IRF8) 
#except PROX1 and TFAP2C, for which the half-site motifs are impossible to generate.  
#Not sure if the PROX1 motif we have in our hand is classified as monomer or homodimer. 

#datas$symbol "FOXD3"   "IRF6"    "MAFF"    "NFATC1"  "ONECUT1" "SOX11"   "TGIF2LX"

#Version2.2

#HT-SELEX
path="../../Data/SELEX-motif-collection/artificial-half-site-motifs/"
length(dir(paste0(path, "orig_pfms/"))) #7


pfms_tab_path=paste0(path,"pfms_tab")
pfms_space_path=paste0(path ,"pfms_space")
pfms_transfac_path=paste0(path,"pfms_transfac")
pfms_scpd=paste0(path,"pfms_scpd") #this can be converted to meme format
pwms_tab_path=paste0(path,"pwms_tab")
pwms_space_path=paste0(path ,"pwms_space")


  
dir.create(pfms_tab_path, recursive=TRUE)
dir.create(pfms_space_path, recursive=TRUE)
dir.create(pfms_transfac_path, recursive=TRUE)
dir.create(pfms_scpd, recursive=TRUE)
dir.create(pwms_tab_path, recursive=TRUE)
dir.create(pwms_space_path, recursive=TRUE)


logos_png_prob=paste0(path, "Logos/png/prob")
logos_png_ic=paste0(path, "Logos/png/ic")
logos_pdf_prob=paste0(path, "Logos/pdf/prob")
logos_pdf_ic=paste0(path, "Logos/pdf/ic")

dir.create(logos_png_prob, recursive=TRUE)
dir.create(logos_png_ic, recursive=TRUE)
dir.create(logos_pdf_prob, recursive=TRUE)
dir.create(logos_pdf_ic, recursive=TRUE)


rev_comp_pfms_tab_path=paste0(path, "PWMs-reverse-complement/pfms_tab")
rev_comp_pfms_space_path=paste0(path, "PWMs-reverse-complement/pfms_space")

dir.create(rev_comp_pfms_tab_path, recursive=TRUE)
dir.create(rev_comp_pfms_space_path, recursive=TRUE)
dir.create(paste0(path, "PWMs-reverse-complement/pfms_transfac"), recursive=TRUE) #Is this needed?

revcomp_logos_png_prob=paste0(path, "PWMs-reverse-complement/Logos/png/prob")
revcomp_logos_pdf_prob=paste0(path, "PWMs-reverse-complement/Logos/pdf/prob")
revcomp_logos_png_ic=paste0(path, "PWMs-reverse-complement/Logos/png/ic")
revcomp_logos_pdf_ic=paste0(path, "PWMs-reverse-complement/Logos/pdf/ic")


dir.create(revcomp_logos_png_prob, recursive=TRUE)
dir.create(revcomp_logos_png_ic, recursive=TRUE)
dir.create(revcomp_logos_pdf_prob, recursive=TRUE)
dir.create(revcomp_logos_pdf_ic, recursive=TRUE)



datas=as_tibble(do.call(rbind,
                       strsplit(
                         dir(paste0(path, "orig_pfms/")),
                         split="_")))


datas$V7=gsub(".pfm", "", datas$V7)
names(datas)=c("symbol","ligand", "batch", "seed", "multinomial", "cycle", "short")
datas$study="artificial-halfsites"
datas$type=NA
datas$organism="Homo_sapiens"
datas$clone=NA
datas$ID=NA

datas$family=NA
datas$experiment="HT-SELEX"
datas$representative=NA
datas$comment=NA

datas$filename=NA
  
datas=datas[, c("ID","symbol",	"clone","family",	"organism","study","experiment",
                "ligand",	"batch", "seed",	"multinomial",
                "cycle","short","representative", "type","comment","filename")]

append=FALSE


datas$IC=NA
datas$length=NA

datas$consensus=""
datas$kld_between_revcomp=NA #Kullback-Leibler divergence between the motif and its reverse complement. Gives idea whether the motif is palindromic.


#Remove these from metadata, empty TFs
remove=c()


for(m in 1:nrow(datas)){
  print(m)
  f=paste0(path,"orig_pfms/" ,dir(paste0(path, "orig_pfms/")))[m]
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
  
    if(datas$experiment[m]=="HT-SELEX"){
      filename=paste0(strsplit(f, "/")[[1]][-c(1,2,3,4,5,6)],collapse="/")  
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
                file=paste0(pfms_tab_path,"/",filename_final) ,sep="\t")
    
  
    write.table(PWM,row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(pfms_space_path,"/",filename_final) ,sep=" ")
    
    
    pcm=as.matrix( read.csv(paste0(pfms_tab_path,"/",filename_final), header=FALSE, sep="\t") )
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
    
    
    pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                profileMatrix=pcm, name=datas$symbol[m],
    )
    
    pfm_rev_comp <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),name=datas$symbol[m],
                                         profileMatrix=pcm_rev_comp
    )
    
    pwm_class <- TFBSTools::toPWM(pfm, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
    
    pwm_class_revcomp <- TFBSTools::toPWM(pfm_rev_comp, type="prob", 
                                          pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
    
    #Compute the KL-divergence between the original and the reverse complement
    
    # Compute KL divergence, it does not matter in which order you give the matrices to the function
    kl_divergence <- as.numeric(compute_kl_divergence(pwm_class@profileMatrix, pwm_class_revcomp@profileMatrix))
    
    datas[m, "kld_between_revcomp"]=kl_divergence
  
    icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
   
    datas$IC[m]=sum(rowSums(icm)) 
    
    #motif length
    datas$length[m]=length(pfm)
    
    motif=universalmotif::read_matrix(file=paste0(pfms_tab_path,"/",filename_final), sep="\t", header=FALSE)
    
    datas[m, "consensus"]=motif@consensus
    
    
    #Write transfac format
    
    transfac=paste0(pfms_transfac_path,"/",filename_final)
    
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
    
    write.table(pwm_class@profileMatrix, file=paste0(pwms_space_path,"/",datas$ID[m], ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
    write.table(pwm_class@profileMatrix, file=paste0(pwms_tab_path,"/",datas$ID[m], ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
    
    
    #Write reverse complements
    
    write.table(pcm_rev_comp, file=paste0(rev_comp_pfms_space_path,"/",datas$ID[m], ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
    write.table(pcm_rev_comp, file=paste0(rev_comp_pfms_tab_path,"/",datas$ID[m], ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
    
    
  
    #Write .scpd format, all into a single file
    PWM=as.matrix(PWM, dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
                append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(pfms_scpd,"/","all.scpd"))
    append=TRUE
    
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
                file=paste0(pfms_scpd,"/","all.scpd"))
    
    datas$filename[m]=paste0(pfms_tab_path,"/",filename_final)
    
    
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

# Merge two data frames
datas <- datas %>%
  left_join(select(TF_family_mappings, "symbol", "Gene Information/ID","DBD" ), by = 'symbol')

datas<- add_column(datas, Lambert2018_families = datas$DBD, .before = 4)


datas <- datas %>% #Save this info only for human monomers
  dplyr::rename("Human_Ensemble_ID"= "Gene Information/ID") %>%
  relocate("Human_Ensemble_ID", .after = "symbol")


column_order=c( "ID", "symbol", "Human_Ensemble_ID","clone","family", "Lambert2018_families", "organism", "study","experiment",          
                "ligand", "batch","seed", "multinomial","cycle","representative",  "short",               
                "type", "comment",  "filename", "IC",  "length","consensus", "kld_between_revcomp")  

datas=datas[,column_order]



write.table(datas, file=paste0(path, "metadata.csv"), row.names = FALSE,sep="\t")


write.table(gsub("../../","",datas$filename), 
            file=paste0(path, "filenames.csv"),
            row.names = FALSE, col.names=FALSE, quote=FALSE)
#save motifnames
write.table(datas$ID, 
            file=paste0(path, "motifnames.csv"),row.names = FALSE, col.names=FALSE, quote=FALSE)





