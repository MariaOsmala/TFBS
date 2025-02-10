library(readr)
library(tidyverse)
library(ggplot2)

k <- as.numeric(commandArgs(trailingOnly = TRUE)) #1-1326

#Create all possible k-mers of length 1-17

bases <- c("A", "C", "G", "T")

# Function to generate k-mers for a given k
generate_kmers <- function(k, bases) {
  print(k)
  if (k == 1) {
    return(bases)
  } else {
    # Use expand.grid to create all combinations of length k
    kmers <- do.call(paste0, expand.grid(rep(list(bases), k)))
    return(kmers)
  }
}

# generate_kmers <- function(k, bases) {
#   print(k)
#   if (k == 1) {
#     return(bases)
#   } else {
#     # Use expand.grid to create all combinations of length k
#     kmers <- apply(permutations(n=length(bases),r=k,v=bases,repeats.allowed=T),1, paste, collapse="")
#     return(kmers)
#   }
# }

# generate_kmers <- function(k, bases) {
#   print(k)
#   if (k == 1) {
#     return(bases)
#   } else {
#     kmers <- apply(unique(t(combn(rep(bases, k), m = k))),1, paste, collapse="")
#     return(kmers)
#   }
# }


# Generate k-mers for lengths 1 to 17 and store in a list
all_kmers <- lapply(k, generate_kmers, bases)


saveRDS(all_kmers, file=paste0("/scratch/project_2006203/TFBS/RProjects/TFBS/RData/",k,"_mers.RDS"))