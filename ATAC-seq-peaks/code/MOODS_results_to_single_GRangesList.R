library(rtracklayer)
library(GenomicRanges)

data_path="/scratch/project_2006203/TFBS/"

#Read motif matches
#dir("../../Results/MOODS_RDS/")

#files=list.files("../../Results/MOODS_RDS/",pattern = "_smaller") #399
#files=list.files("../../Results/MOODS_Vierstra_processed/MOODS_Vierstra_RDS/") #399


#files=list.files(paste0(data_path, "Results/MOODS_Teemu_processed/MOODS_RDS"))
#files=list.files(paste0(data_path, "Results/MOODS_Teemu_4_processed/MOODS_RDS"))
#files=list.files(paste0(data_path, "Results/MOODS_Teemu_3_processed/MOODS_RDS"))
#files=list.files(paste0(data_path, "Results/MOODS_Teemu_2_processed/MOODS_RDS"))
#files=list.files(paste0(data_path, "Results/MOODS_Teemu_rest_processed/MOODS_RDS"))
#files=list.files(paste0(data_path, "Results/MOODS_Mouse_processed/MOODS_RDS"))
files=list.files(paste0(data_path, "Results/MOODS_Mouse_processed_4/MOODS_RDS"))

files=files[grep("top", files)]

GR_list=GRangesList()

for(file in files){
  print(file)
  #GR_list=c(GR_list, readRDS(paste0("../../Results/MOODS_RDS/",file)) )
  #GR_list=c(GR_list, readRDS(paste0("../../Results/MOODS_Vierstra_processed/MOODS_Vierstra_RDS/",file)) )
  
  #GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_Teemu_processed/MOODS_RDS/",file)) )
  #GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_Teemu_4_processed/MOODS_RDS/",file)) )
  #GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_Teemu_3_processed/MOODS_RDS/",file)) )
  #GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_Teemu_2_processed/MOODS_RDS/",file)) )
  #GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_Teemu_rest_processed/MOODS_RDS/",file)) )
  #GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_Mouse_processed/MOODS_RDS/",file)) )
  GR_list=c(GR_list, readRDS(paste0(data_path,"Results/MOODS_Mouse_processed_4/MOODS_RDS/",file)) )
  
}

#saveRDS(GR_list,  file = paste0( "../RData/all_motif_matches_sbatch.Rds")) #
#saveRDS(GR_list,  file = paste0( "../RData/all_motif_matches_Vierstra.Rds")) #
#saveRDS(GR_list,  file = paste0( "../RData/top_motif_matches_Teemu.Rds")) #
#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_Teemu.Rds")) #
#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_4_Teemu.Rds")) #
#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_3_Teemu.Rds")) #
#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_2_Teemu.Rds")) #
#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_Teemu_rest.Rds")) #
#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_Mouse.Rds")) #
saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_Mouse_4.Rds")) #
