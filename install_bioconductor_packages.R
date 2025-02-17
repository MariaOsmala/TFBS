# Read the YAML file
#setwd("~/projects/TFBS")
library(yaml)
config <- yaml.load_file("Experiments/config.yaml", readLines.warn=FALSE)

# Print the loaded configuration
# print(config)

.libPaths(config$R_library_dir) #Or add this to ~/.Renviron

# Check if BiocManager is installed, and install it if missing

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::version()
#‘3.20’

# List of Bioconductor packages to install

# List of Bioconductor packages to install
bioconductor_packages <- c(
  "BiocGenerics",
  "rtracklayer",
  "ensembldb",
  "org.Hs.eg.db",
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Hsapiens.UCSC.hg38",
  "GenomicDistributions",
  "ChIPseeker",
  "universalmotif",
  # "GeneTonic",
  "genomation",
  "motifmatchr",
  "motifStack",
  "GenomicFeatures",
  "Rsamtools",
  "monaLisa",
  "GO.db",
  "BiocFileCache",
  "annotate",
  "TFBSTools",
  "TFMPvalue",
  "seqPattern",
  "KEGGREST",
  "impute",
  "BiocVersion",
  "DirichletMultinomial",
  "BiocParallel",
  "CNEr"
)

Sys.getenv("CC")
Sys.getenv("CXX")

# Install missing Bioconductor packages
for (pkg in bioconductor_packages) {
  #pkg=bioconductor_packages[1]
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)  # ask=FALSE prevents prompts
  } else {
    message(paste("Package", pkg, "is already installed."))
  }
}

#- r-webgestaltr

# Install commonly used CRAN packages
#cran_packages <- c(
#  "r-webgestaltr"
#)

#install.packages(setdiff(cran_packages, rownames(installed.packages())))
devtools::install_github("bzhanglab/WebGestaltR") #Needs Rust

# Save session info to a text file
sink(paste0(config$base_dir, "session-info-r-4.4.2.txt"))
sessionInfo()
sink()



# Confirm installation
message("All Bioconductor and CRAN packages installed successfully!")

