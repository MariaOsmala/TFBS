
library("rtracklayer")
library("dbplyr")
library("dplyr")

library("ggplot2")

setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")



data_path="/scratch/project_2006203/TFBS/"

representatives=read.table("../../PWMs_final_version2.2/metadata_representatives.tsv",sep="\t", header=TRUE) #3933

rep_motifs=representatives$ID


#grep("CTCF",rep_motifs) CTCF is in representative motif

length(unique(rep_motifs))
#3933

#Is the order here same as in the representatives, this is 56 GB, seems so
#representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_human_final_version2.2.Rds")

arrays=1:394

#motifs_all=names(representative_motif_matches)
motifs_all=representatives$ID

cCREs_list<- readRDS(file = paste0(data_path, "ATAC-seq-peaks/RData/cCREs_list.Rds") ) #
length(cCREs_list)

p_matrix=matrix(nrow=length(motifs_all), ncol=length(cCREs_list),
                dimnames = list( motifs_all,
                                 names(cCREs_list) )) #3933 x 111

e_matrix=matrix(nrow=length(motifs_all), ncol=length(cCREs_list),
                dimnames = list(motifs_all,
                                names(cCREs_list) )) #3933 x 111



for(array in arrays){
  print(array)
  #array=1
  #start=(array-1)*10+1
  #end=(array)*10
  
  #if(end>3933){
  #  end=3933
  #}
  
  motifs=motifs_all[start:end]
  
  p_matrix_parts=readRDS( paste0(data_path, "ATAC-seq-peaks/RData/final_p_matrix_human_cell_type_restricted_all_motifs_version2.2_array_",array,".Rds"))
  e_matrix_parts=readRDS( paste0( data_path, "ATAC-seq-peaks/RData/final_e_matrix_human_cell_type_restricted_all_motifs_version2.2_array_",array,".Rds")) 
  
  p_matrix[rownames(p_matrix_parts),]=p_matrix_parts
  e_matrix[rownames(e_matrix_parts),]=e_matrix_parts
  
}


saveRDS(p_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/final_p_matrix_human_cell_type_restricted_all_motifs_version2.2_fromParts.Rds")) #
saveRDS(e_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/final_e_matrix_human_cell_type_restricted_all_motifs_version2.2_fromParts.Rds")) #





