library(yaml)
config <- yaml.load_file("Experiments/config.yaml", readLines.warn=FALSE)

# Print the loaded configuration
# print(config)

.libPaths(config$R_library_dir)


setwd(paste0(config$base_dir, "code/SELEX-motif-collection"))

#Draws logos using ggseqlogo package

source("extract_motifs_from_excel_Jolma2013.R") 
source("extract_motifs_Morgunova_2015.R") 
source("extract_motifs_from_excel_Jolma2015.R") 
source("extract_motifs_from_excel_Nitta2015.R") 
source("extract_motifs_from_excel_Yin2017.R") 
source("extract_motifs_from_excel_Xie2025.R") 


#Creates Data/SELEX-motif-collection/metadata.tsv, note that does not contain the representativenes info
source("collect_all_motif_metadata.R")  

#Artificial half-sites, not included in SELEX collection but needed in some analysis
source("extract_motifs_from_excel_artificial_halfsites.R")


#run this first
#~/projects/TFBS/code/SELEX-Dominating-Set-Analysis/process_DSA_results.R
#then run
#source("add_match_numbers_and_match_scores.R")

#representatives separately