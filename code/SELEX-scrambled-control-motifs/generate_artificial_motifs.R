
library(readr)
rm(list=ls())

source("artificial_motif_functions.R")
source("../SELEX-motif-collection/paths_to_folders.R")


metadata <- read_delim("../../Data/SELEX-motif-collection/metadata.csv",
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)


artificial_all<-list()



for(i in 1:nrow(metadata) ){
  #i=1
  print(i)
  pcm=as.matrix( read.table(paste0( metadata$filename[i]), header=FALSE) )
  
  dimnames(pcm)=list(c("A", "C", "G", "T"))
  
  min_slice_width = floor( ncol(pcm) / 3 )
  
  artificial<-list()
  for( j in min_slice_width:(ncol(pcm) - min_slice_width) ){
    #j=5
    # Get matrix versions where slices are concatenated in all 3 relative orientations
   
    artificial[paste0(metadata$ID[i], "_",j,c("_HT2", "_HH", "_TT")) ]=get_other_pfm_slice_pair_orientations(pcm, j)
    
  }
  
  for(motif in names(artificial)){
    print(motif)
    write.table(artificial[[motif]],row.names = FALSE, col.names=FALSE, quote=FALSE,file=paste0(scrambled_pfms_tab_path,"/", motif, ".pfm"), sep="\t")
    write.table(artificial[[motif]],row.names = FALSE, col.names=FALSE, quote=FALSE,file=paste0(scrambled_pfms_space_path,"/", motif, ".pfm"), sep=" ")
    
    
    transfac=universalmotif::read_matrix(file=paste0(scrambled_pfms_tab_path,"/", motif, ".pfm"), sep="\t", header=FALSE)
    transfac@name=motif
    

    transfac_file=paste0(scrambled_pfms_tab_path,"/", motif,".pfm")
    universalmotif::write_transfac(transfac, file=transfac_file, overwrite = TRUE, append = FALSE)
    
    
  }
  
  artificial_all[[metadata$ID[i]]]=artificial

}


saveRDS(artificial_all, file="../../RData/scrambled_control_motifs.Rds")






