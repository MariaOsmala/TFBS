
library("readr")
library("dplyr")

df_motif_info=read_delim("../../Data/SELEX-motif-collection/metadata_representative_matchNumber_threshold_phyloP.tsv", "\t")


df_motif_info %>% pull(new_representative) %>% table()

#NO  YES 
#2701 1232 

df_motif_info$representative=NULL

#rename column 

colnames(df_motif_info)[colnames(df_motif_info)=="new_representative"]="representative"

df_motif_info %>% pull(representative) %>% table()


table(df_motif_info$Methyl.SELEX.Motif.Category, useNA = "always")

# Inconclusive       MethylPlus Multiple effects             <NA> 
#   31              245               21             3636 
#31      +        245      +         21 = 297

write_delim(df_motif_info, "../../Data/SELEX-motif-collection/metadata_final.tsv", "\t", col_names = TRUE)

df_motif_info=df_motif_info %>% filter(representative=="YES")


output_folder="../../Data/forSharing/representatives/"

#create folders 

dir.create(paste0(output_folder,"pfms_space"), showWarnings = FALSE)
dir.create(paste0(output_folder,"pfms_tab"), showWarnings = FALSE)
dir.create(paste0(output_folder,"pfms_transfac"), showWarnings = FALSE)

#copy files to another folder

for(i in 1:nrow(df_motif_info)){
  
 file.copy(from = df_motif_info$filename[i], to = paste0(output_folder,"pfms_tab/",strsplit(df_motif_info$filename[i], "/")[[1]][6] ))
  
 file.copy(from = gsub("_tab", "_space", df_motif_info$filename[i]), to = paste0(output_folder,"pfms_space/",strsplit(df_motif_info$filename[i], "/")[[1]][6] ))
  
 #transfac format has some number in the filename?
 file.copy(from = gsub("_tab", "_transfac", df_motif_info$filename[i]), to = paste0(output_folder,"pfms_transfac/",strsplit(df_motif_info$filename[i], "/")[[1]][6] ))
  
}


write_delim(df_motif_info, paste0(output_folder, "metadata.tsv"), "\t", col_names = TRUE)
