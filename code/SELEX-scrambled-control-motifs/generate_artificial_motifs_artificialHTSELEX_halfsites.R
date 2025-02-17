rm(list=ls())

source("~/projects/TFBS/RProjects/TFBS/code/artificial_motif_functions.R")

library(readr)

setwd("/Users/osmalama/projects/TFBS/RProjects/TFBS")


# Generate the following folders  and save the motifs in them


if (!dir.exists("~/projects/TFBS/PWMS_final_version2.2/fromYimeng/half-site-motifs-20240426/artificial_motifs/")) {
  dir.create("~/projects/TFBS/PWMS_final_version2.2/fromYimeng/half-site-motifs-20240426/artificial_motifs/")
}

if (!dir.exists("~/projects/TFBS/PWMS_final_version2.2/fromYimeng/half-site-motifs-20240426/artificial_motifs_space/")) {
  dir.create("~/projects/TFBS/PWMS_final_version2.2/fromYimeng/half-site-motifs-20240426/artificial_motifs_space/")
}

if (!dir.exists("~/projects/TFBS/PWMS_final_version2.2/fromYimeng/half-site-motifs-20240426/artificial_motifs_transfac/")) {
  dir.create("~/projects/TFBS/PWMS_final_version2.2/fromYimeng/half-site-motifs-20240426/artificial_motifs_transfac/")
}

#Need to also generate individual motifs of the half-sites, but need to check 
#that these do not correspond to any true motifs?


#metadata_union <- read_delim("~/projects/TFBS/PWMs_final_union/metadata_4003_motifs.csv", 
#                           delim = "\t", escape_double = FALSE, 
#                           trim_ws = TRUE)


metadata <- read_delim("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/metadata.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)



artificial_all<-list()

# metadata=metadata[new_ind,]

for(i in 1:nrow(metadata) ){
  #i=1
  print(i)
  pcm=as.matrix( read.table(paste0("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/pwms/", metadata$ID[i], ".pfm"), header=FALSE) )
  
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
    write.table(artificial[[motif]],row.names = FALSE, col.names=FALSE, quote=FALSE,file=paste0("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/artificial_motifs/", motif, ".pfm"), sep="\t")
    write.table(artificial[[motif]],row.names = FALSE, col.names=FALSE, quote=FALSE,file=paste0("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/artificial_motifs_space/", motif, ".pfm"), sep=" ")
    
    
    transfac=universalmotif::read_matrix(file=paste0("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/artificial_motifs/", motif, ".pfm"), sep="\t", header=FALSE)
    transfac@name=motif
    

    transfac_file=paste0("../../PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/artificial_motifs_transfac/", motif,".pfm")
    universalmotif::write_transfac(transfac, file=transfac_file, overwrite = TRUE, append = FALSE)
    
    
  }
  
  artificial_all[[metadata$ID[i]]]=artificial
  
  #Write all motifs to file:
  

}

#saveRDS(artificial_all, file="RData/artificial_all.Rds")
saveRDS(artificial_all, file="RData/artificial_artificialHTSelex_halfsites.Rds")





