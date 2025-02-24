# Read the YAML file
#setwd("~/projects/TFBS")
.libPaths("/Users/osmalama/.local/share/mamba/envs/TFBS-R4.2.1/lib/R/library")
library(yaml)
config <- yaml.load_file("Experiments/config.yaml", readLines.warn=FALSE)

# Print the loaded configuration
# print(config)

.libPaths(config$R4.2.1_library_dir) #Or add this to ~/.Renviron

# Check if BiocManager is installed, and install it if missing

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::version() # or 3.16 if r-4.2.1


#If using R-4.2.1, use version = '3.16'

Sys.getenv("CC")
Sys.getenv("CXX")

Sys.setenv("CC"="/Users/osmalama/.local/share/mamba/envs/TFBS-R4.2.1/bin/clang")
Sys.setenv("CXX"="/Users/osmalama/.local/share/mamba/envs/TFBS-R4.2.1/bin/clang++")


# For R-4.2.1
bioconductor_packages <- c(
  "AnnotationDbi", #OK
  "GenomicRanges", #OK
  "IRanges", #OK
  "S4Vectors", #OK
  "BiocGenerics", #OK
  "GenomeInfoDb", #OK
  "rtracklayer",
  "ensembldb",
  "org.Hs.eg.db",
  "biomaRt",
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Hsapiens.UCSC.hg38",
  "GenomicDistributions",
  "ChIPseeker",
  "Biostrings",
  "AnnotationHub",
  "universalmotif",
  "GeneTonic",
  "genomation",
  "BiocParallel",
  "ComplexHeatmap",
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
    BiocManager::install(pkg, ask = TRUE)  # ask=FALSE prevents prompts
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
sink(paste0(config$base_dir, "session-info-r-4.2.1.txt"))
sessionInfo()
sink()



# Confirm installation
message("All Bioconductor and CRAN packages installed successfully!")







