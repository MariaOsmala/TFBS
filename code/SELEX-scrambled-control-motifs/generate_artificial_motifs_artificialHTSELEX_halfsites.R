rm(list=ls())

source("artificial_motif_functions.R")

library(readr)

scrambled_pfms_tab_path="../../Data/SELEX-motif-collection/artificial-half-site-motifs/scrambled_pfms_tab"
scrambled_pfms_space_path="../../Data/SELEX-motif-collection/artificial-half-site-motifs/scrambled_pfms_space"
scrambled_pfms_transfac_path="../../Data/SELEX-motif-collection/artificial-half-site-motifs/scrambled_pfms_tranfac"
scrambled_pfms_scpd="../../Data/SELEX-motif-collection/artificial-half-site-motifs/scrambled_pfms_scpd" #this can be converted to meme format

scrambled_pwms_tab_path="../../Data/SELEX-motif-collection/artificial-half-site-motifs/scrambled_pwms_tab"
scrambled_pwms_space_path="../../Data/SELEX-motif-collection/artificial-half-site-motifs/scrambled_pwms_space"


dir.create(scrambled_pfms_tab_path, recursive=TRUE)
dir.create(scrambled_pfms_space_path, recursive=TRUE)
dir.create(scrambled_pfms_transfac_path, recursive=TRUE)
dir.create(scrambled_pfms_scpd, recursive=TRUE)
dir.create(scrambled_pwms_tab_path, recursive=TRUE)
dir.create(scrambled_pwms_space_path, recursive=TRUE)


metadata <- read_delim("../../Data/SELEX-motif-collection/artificial-half-site-motifs/metadata.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)



artificial_all<-list()



for(i in 1:nrow(metadata) ){
  #i=1
  print(i)
  pcm=as.matrix( read.table(paste0("../../Data/SELEX-motif-collection/artificial-half-site-motifs/pfms_tab/", metadata$ID[i], ".pfm"), header=FALSE) )
  
  dimnames(pcm)=list(c("A", "C", "G", "T"))
  
  min_slice_width = floor( ncol(pcm) / 3 )
  
  artificial<-list()
  for( j in min_slice_width:(ncol(pcm) - min_slice_width) ){
    artificial[paste0(metadata$ID[i], "_",j,c("_HT2", "_HH", "_TT")) ]=get_other_pfm_slice_pair_orientations(pcm, j)
    
  }
  
  for(motif in names(artificial)){
    print(motif)
    write.table(artificial[[motif]],row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(scrambled_pfms_tab_path,"/", motif, ".pfm"), sep="\t")
    write.table(artificial[[motif]],row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(scrambled_pfms_space_path,"/", motif, ".pfm"), sep=" ")
    
    
    transfac=universalmotif::read_matrix(file=paste0(scrambled_pfms_tab_path,"/", motif, ".pfm"), sep="\t", header=FALSE)
    transfac@name=motif
    

    transfac_file=paste0(scrambled_pfms_transfac_path,"/", motif,".pfm")
    universalmotif::write_transfac(transfac, file=transfac_file, overwrite = TRUE, append = FALSE)
    
    
  }
  
  artificial_all[[metadata$ID[i]]]=artificial
  

}


saveRDS(artificial_all, file="../../RData/scrambled_control_motifs_artificialHTSelex_halfsites.Rds")





