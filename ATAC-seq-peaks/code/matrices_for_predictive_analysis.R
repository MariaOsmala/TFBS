library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")
#library(ggbreak)
#install.packages("ggbreak")
#library("GRanges")
library("readr")

#The logistic regression model predicts whether a cCRE i is specific to certain cell type (Class 1) or not (0). 
#The model is trained for all 111 cell types.

#In the simplest model, the covariates are presence/absence (0 or 1)  or counts (0,1,2,...)  
#of motif matches (individual representative PWMS ~1000 features). 
#Or use the binding affinities of motif matches? Positional logistic regression classifier? 
#Use R-packages caret and glmnet. 


#Database stuff
setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")

#source("../code/FishersExact.R")

data_path="/scratch/project_2006203/TFBS/"


#Motif metadata with ICs and length
#representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294 version 1
representatives=read.table("../../PWMs_final_version2.2/metadata_representatives.tsv",sep="\t", header=TRUE) #3933
rep_motifs=representatives$ID[which(representatives$new_representative=="YES")]

length(unique(rep_motifs))
#1031

representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_human_final_version2.2.Rds")


motifs=names(representative_motif_matches)


cCREs_list<- readRDS(file = paste0(data_path, "ATAC-seq-peaks/RData/cCREs_list.Rds") ) #
length(cCREs_list)

Adult_Celltypes_mapping <- read_delim("/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

#typo
Adult_Celltypes_mapping$`Cell type`[grep("Alverolar Type 2,Immune", Adult_Celltypes_mapping$`Cell type`)]="Alveolar Type 2,Immune"

names(cCREs_list) %in% Adult_Celltypes_mapping$celltype

Adult_Celltypes_mapping$celltype %in% names(cCREs_list) 

#test=data.frame( names(cCREs_list), Adult_Celltypes_mapping$`Cell type`[ match(names(cCREs_list), Adult_Celltypes_mapping$celltype) ])

names(cCREs_list)=Adult_Celltypes_mapping$`Cell type`[ match(names(cCREs_list), Adult_Celltypes_mapping$celltype) ]

sum(sapply(cCREs_list, length)) #1091805

all_cCREs=unlist(cCREs_list) #which are unique 1091805, 
#length(all_cCREs) #1091805, all are of width 400 bp

#table(width(all_cCREs))

#consider only unique
all_cCREs=unique(all_cCREs) #435142


N=length(all_cCREs) #all possible cCREs #435142

presence_matrix=matrix(0,nrow=N,ncol=length(motifs), 
                       dimnames=list(NULL, motifs)) #N x 1000

count_matrix=matrix(0,nrow=N,ncol=length(motifs), 
                    dimnames=list(NULL, motifs)) #N x 1000

#How to take the scores of multiple matches into account
max_score_matrix=matrix(0,nrow=N,ncol=length(motifs), 
                        dimnames=list(NULL, motifs)) #N x 1000



for(hits in motifs){
  print(hits)
  #hits=motifs[1]
  
  fol=findOverlaps(representative_motif_matches[[hits]], all_cCREs, type="within", ignore.strand=TRUE)
  presence_matrix[unique(fol@to),hits]=1
  count_matrix[as.numeric(names(table(fol@to))), hits]=as.vector(table(fol@to))
  
  indices=lapply(unique(fol@to), function(x, y) which(y==x), fol@to)
  matches=lapply(indices, function (x,y,z) y[z[x]] , y=representative_motif_matches[[hits]], z=fol@from)
  max_score_matrix[unique(fol@to), hits]=sapply(matches, function(x) max(x$score))
  
  
  
}

#saveRDS(presence_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/presence_matrix.Rds")) #version 1
#saveRDS(count_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/count_matrix.Rds")) #version 1
#saveRDS(max_score_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/max_score_matrix.Rds")) #version 1


saveRDS(presence_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/presence_matrix_version2.2.Rds")) #version 2.2
saveRDS(count_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/count_matrix_version2.2.Rds")) #version 2.2
saveRDS(max_score_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/max_score_matrix_version2.2.Rds")) #version 2.2





