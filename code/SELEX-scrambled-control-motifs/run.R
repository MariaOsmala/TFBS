library(yaml)
config <- yaml.load_file("Experiments/config.yaml", readLines.warn=FALSE)

# Print the loaded configuration
# print(config)

.libPaths(config$R_library_dir)


setwd(paste0(config$base_dir, "code/SELEX-scrambled-control-motifs/"))


source("generate_artificial_motifs.R")

source("generate_artificial_motifs_artificialHTSELEX_halfsites.R")
