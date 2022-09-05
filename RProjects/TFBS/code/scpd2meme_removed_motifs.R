# ./PWMs/Jolma2015/Homo_sapiens_all.meme
# Error: Motif position 11 summed to zero.
# Error: Motif position 11 summed to zero.
# Error: Motif position 11 summed to zero.
# Error: Motif position 11 summed to zero.
# Error: Motif position 9 summed to zero.
# Error: Motif position 11 summed to zero.
# Error: Motif position 10 summed to zero.
# Error: Motif position 11 summed to zero.
# Error: Motif position 11 summed to zero.
# Error: Motif position 11 summed to zero.
# Converted 652 motifs.
# Skipped 0 motifs.
# 10 conversion errors.

metadata <- readRDS("/home/local/osmalama/projects/TFBS/RProjects/TFBS/data/Jolma2015.Rds")

motifnames <- apply(metadata[,-which(colnames(metadata)%in% c("clone", "family", "organism", "study","comment", "short", "type", "filename"))], 1,paste0,collapse="_")
length(motifnames)
library(tidyverse)
library(readxl)
source("code/read_excel_tables.R")

meme<-read.table("../../PWMs/Jolma2015/Homo_sapiens_all.meme", fill=TRUE)
length(meme[which(meme$V1=="MOTIF"),2])

length(which(motifnames %in% meme[which(meme$V1=="MOTIF"),2]))
missing=motifnames[-which(motifnames %in% meme[which(meme$V1=="MOTIF"),2])]
#10 motifs for which the last columns are zero, these are not included in meme
#

# ./PWMs/fromYimeng/Homo_sapiens_all.meme
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Converted 695 motifs.
# Skipped 0 motifs.
# 10 conversion errors.

metadata <- readRDS("/home/local/osmalama/projects/TFBS/RProjects/TFBS/data/fromYimeng.Rds")

motifnames <- apply(metadata[,-which(colnames(metadata)%in% c("clone", "family", "organism", "study","experiment","comment","representative",  "filename"))], 1,paste0,collapse="_")
length(motifnames)
library(tidyverse)

meme<-read.table("../../PWMs/fromYimeng/Homo_sapiens_all.meme", fill=TRUE)
length(meme[which(meme$V1=="MOTIF"),2])

length(which(motifnames %in% meme[which(meme$V1=="MOTIF"),2]))
missing=motifnames[-which(motifnames %in% meme[which(meme$V1=="MOTIF"),2])]
#these are all zero matrices


# 
# ./PWMs/Yin2017/Homo_sapiens_all.meme
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Empty matrix.
# Error: Empty matrix.
# Error: Empty matrix.
# Error: Empty matrix.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Error: Expected count matrix but got probability matrix. Can't continue without site count.
# Converted 1787 motifs.
# Skipped 0 motifs.
# 7 conversion errors.
metadata <- readRDS("/home/local/osmalama/projects/TFBS/RProjects/TFBS/data/Yin2017.Rds")

motifnames <- apply(metadata[,-which(colnames(metadata)%in% c("clone", "family", "organism", "study","comment", "short", "type", "filename"))], 1,paste0,collapse="_")
length(motifnames)
library(tidyverse)

meme<-read.table("../../PWMs/Yin2017/Homo_sapiens_all.meme", fill=TRUE)
length(meme[which(meme$V1=="MOTIF"),2])

length(which(motifnames %in% meme[which(meme$V1=="MOTIF"),2]))
missing=motifnames[-which(motifnames %in% meme[which(meme$V1=="MOTIF"),2])]
#zero or empty matrices
