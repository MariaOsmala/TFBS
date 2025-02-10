library(rtracklayer)
library(GenomicRanges)

data_path="/scratch/project_2006203/TFBS/"

#Read motif matches


#files=list.files(paste0(data_path, "Results/MOODS_human_final_processed/MOODS_RDS"))
files=list.files(paste0(data_path, "Results/MOODS_human_final_version2.2_correct_processed/MOODS_RDS"))
#files=list.files(paste0(data_path, "Results/MOODS_human_final_artificialHTSelex_processed/MOODS_RDS"))
#files=list.files(paste0(data_path, "Results/MOODS_mouse_mm39_final_processed/MOODS_RDS"))


files=files[grep("top", files)]
GR_list=GRangesList()

for(file in files){
  print(file)
  
  
  #GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_human_final_processed/MOODS_RDS/",file)) )
  GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_human_final_version2.2_correct_processed/MOODS_RDS/",file)) )
  #GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_human_final_artificialHTSelex_processed/MOODS_RDS/",file)) )
  #GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_mouse_mm39_final_processed/MOODS_RDS/",file)) )
}



#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_human_final.Rds")) #
saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_human_final_version2.2.Rds")) #
#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_human_artificialHTSelex_version2.2.Rds")) #

#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_mouse_final.Rds")) #
