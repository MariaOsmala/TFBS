library(rtracklayer)
library(GenomicRanges)

data_path="/scratch/project_2006203/TFBS/"

#Read motif matches

args <- commandArgs(trailingOnly = TRUE)
MOODS_path=args[1]
results_path=args[2]


#MOODS_path="Results/MOODS_human_final_version2.2_correct_processed/MOODS_RDS"
#results_path="/ATAC-seq-peaks/RData/top_motif_matches_human_final_version2.2.Rds"


files=list.files(paste0(data_path, MOODS_path))

files=files[grep("top", files)]
GR_list=GRangesList()

for(file in files){
  print(file)
  GR_list=c(GR_list, readRDS(paste0(data_path,MOODS_path,file)) )
  
}


saveRDS(GR_list,  file = paste0(data_path, results_path)) #
