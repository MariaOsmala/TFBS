library(rtracklayer)
library(GenomicRanges)

data_path="/scratch/project_2006203/TFBS/"

#Read motif matches

results_path=paste0(data_path, "Results/MOODS_human_final_version2.2_processed/MOODS_bigbed")

#files=list.files(paste0(data_path, "Results/MOODS_human_final_processed/MOODS_RDS"))
#files=list.files(paste0(data_path, "Results/MOODS_human_final_union_processed/MOODS_RDS"))
#files=list.files(paste0(data_path, "Results/MOODS_mouse_mm39_final_processed/MOODS_RDS"))
files=list.files(results_path)


#files=files[grep("top", files)]
GR_list=GRangesList()

for(file in files){
  #file=files[1]
  print(file)
  
  motif=strsplit(file, "_top.bed")[[1]]
  
 
  tmp_GRanges_top=import.bed( paste0(results_path, "/", file))
  
  GR_list[[motif]]=tmp_GRanges_top
  
}



#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_human_final.Rds")) #
#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_human_final_union.Rds")) #
saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_human_final_version2.2.Rds")) #

#saveRDS(GR_list,  file = paste0(data_path, "/ATAC-seq-peaks/RData/top_motif_matches_mouse_final.Rds")) #
