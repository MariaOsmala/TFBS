library(yaml)
config <- yaml.load_file("Experiments/config.yaml", readLines.warn=FALSE)

# Print the loaded configuration
# print(config)

.libPaths(config$R_library_dir)


setwd(paste0(config$base_dir, "code/SELEX-motif-collection"))


source("extract_motifs_from_excel_Jolma2013.R")
source("extract_motifs_Morgunova_2015.R")
source("extract_motifs_from_excel_Jolma2015.R")
source("extract_motifs_from_excel_Nitta2015.R")
source("extract_motifs_from_excel_Yin2017.R")
source("extract_motifs_from_excel_Xie2025.R")

source("collect_all_motif_metadata.R")  


source("extract_motifs_from_excel_artificial_halfsites.R")

