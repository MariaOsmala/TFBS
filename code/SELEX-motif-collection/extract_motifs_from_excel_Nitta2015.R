

#############Read excel tables################

library(tidyverse)
library(readxl)
library(universalmotif)
library(ggplot2)
rm(list=ls())
source("read_excel_tables.R")
source("paths_to_folders.R", echo=FALSE)
source("compute_kl_divergence.R", echo=FALSE)
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


all_table <- read_excel("../../Data/SELEX-motif-collection/elife-04837-supp1-v1.xlsx", 
                        sheet="E. PWMs", col_names=FALSE, skip=2094)

 
PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(seq(1,nrow(PWMs),5), function(i) PWMs[i,]))

PWMs_metadata=cbind(PWMs_metadata[,1], do.call(rbind, sapply(PWMs_metadata[,2], strsplit, split="_")) )

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

PWMs_metadata$ID=""
PWMs_metadata$IC=NA
PWMs_metadata$length=NA
PWMs_metadata$consensus=""
PWMs_metadata$kld_between_revcomp=NA #Kullback-Leibler divergence between the motif and its reverse complement. Gives idea whether the motif is palindromic.




PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

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
    
    filename=paste0(pfms_tab_path,"/", 
                    paste0(PWMs_metadata[m,
                                         -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", 
                                                                              "representative", "short", "type", "filename","ID", "IC","length", "consensus", "kld_between_revcomp"))], collapse="_"),
                    ".pfm")
    ID=paste0(PWMs_metadata[m,
                            -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative", "short", "type", "filename",
                                                                 "ID", "IC", "length", "consensus", "kld_between_revcomp" ))], collapse="_")
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
    PWMs_metadata[m, "IC"]=sum(rowSums(icm)) #6.77892 # universalmotif computes this wrong?
    
    #motif length
    PWMs_metadata[m, "length"]=length(pfm)
    motif=universalmotif::read_matrix(file=filename, sep="\t", header=FALSE)
    motif@name=ID
    
    PWMs_metadata[m, "consensus"]=motif@consensus
  
    transfac=paste0(pfms_transfac_path, "/", ID,".pfm")
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
  
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0(pfms_space_path,"/", ID, ".pfm"), sep=" ")
  
    PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
                append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(pfms_scpd,"/", "Nitta2015",".scpd"))
    append=TRUE
    
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
                file=paste0(pfms_scpd, "/","Nitta2015",".scpd"))
  
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

TF_family_mappings <- read_csv("../../Data/SELEX-motif-collection/HumanTFs-ccbr-TableS1.csv")


# Merge two data frames
PWMs_metadata <- PWMs_metadata %>%
  left_join(select(TF_family_mappings, "symbol", "Gene Information/ID","DBD" ), by = 'symbol')

PWMs_metadata <- add_column(PWMs_metadata, Lambert2018_families = PWMs_metadata$DBD, .before = 4)

PWMs_metadata$DBD=NULL




PWMs_metadata<- PWMs_metadata %>% #Save this info only for human monomers
  dplyr::rename("Human_Ensemble_ID"="Gene Information/ID") %>%
  relocate("Human_Ensemble_ID", .after = "symbol")

PWMs_metadata<- PWMs_metadata %>% #Save this info only for human monomers
  relocate("ID", .before = "symbol")

write.table(PWMs_metadata, file="../../Data/SELEX-motif-collection/Nitta2015_metadata.csv", row.names = FALSE, sep="\t")




