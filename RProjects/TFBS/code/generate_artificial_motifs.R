rm(list=ls())

source("~/projects/TFBS/RProjects/TFBS/code/artificial_motif_functions.R")

library(readr)
#metadata <- read_delim("~/projects/TFBS/PWMs_final/metadata.csv", 
#                           delim = "\t", escape_double = FALSE, 
#                           trim_ws = TRUE)

#Previous artificial motifs are in 
#~/projects/TFBS/artificial_motifs/
#~/projects/TFBS/artificial_motifs_space/
#~/projects/TFBS/artificial_motifs_transfac/

# Generate the following folders  and save the motifs in them


if (!dir.exists("~/projects/TFBS/PWMS_final_version2.2/artificial_motifs/")) {
  dir.create("~/projects/TFBS/PWMS_final_version2.2/artificial_motifs/")
}

if (!dir.exists("~/projects/TFBS/PWMS_final_version2.2/artificial_motifs_space/")) {
  dir.create("~/projects/TFBS/PWMS_final_version2.2/artificial_motifs_space/")
}

if (!dir.exists("~/projects/TFBS/PWMS_final_version2.2/artificial_motifs_transfac/")) {
  dir.create("~/projects/TFBS/PWMS_final_version2.2/artificial_motifs_transfac/")
}

#Need to also generate individual motifs of the half-sites, but need to check 
#that these do not correspond to any true motifs?


#metadata_union <- read_delim("~/projects/TFBS/PWMs_final_union/metadata_4003_motifs.csv", 
#                           delim = "\t", escape_double = FALSE, 
#                           trim_ws = TRUE)


metadata <- read_delim("~/projects/TFBS/PWMs_final_version2.2/metadata_3993_motifs.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)


#Move those artificial motifs that correspond to union that are also in the final version2.2

# earlier_ind=which(metadata$ID %in% metadata_union$ID) #3867
# 
# new_ind=which(!(metadata$ID %in% metadata_union$ID)) #66
# 
# which(metadata$ID=="ZBTB2_HT-SELEX_TTAAAT40NAAT_KT_NTTTMCGGTWAN_1_4")
# 
# for(ei in earlier_ind){
#   #ei=earlier_ind[1]
#   files=grep(metadata$ID[ei], list.files("~/projects/TFBS/PWMS_final_union/artificial_motifs/"), value=TRUE)
#   file.copy(paste0("~/projects/TFBS/PWMS_final_union/artificial_motifs/", files), "~/projects/TFBS/PWMS_final_version2.2/artificial_motifs/", overwrite = TRUE )
#   
#   files=grep(metadata$ID[ei], list.files("~/projects/TFBS/PWMS_final_union/artificial_motifs_space/"), value=TRUE)
#   file.copy(paste0("~/projects/TFBS/PWMS_final_union/artificial_motifs_space/", files), "~/projects/TFBS/PWMS_final_version2.2/artificial_motifs_space/", overwrite = TRUE )
#   
#   files=grep(metadata$ID[ei], list.files("~/projects/TFBS/PWMS_final_union/artificial_motifs_transfac/"), value=TRUE)
#   file.copy(paste0("~/projects/TFBS/PWMS_final_union/artificial_motifs_transfac/", files), "~/projects/TFBS/PWMS_final_version2.2/artificial_motifs_transfac/", overwrite = TRUE )
# }

artificial_all<-list()

# metadata=metadata[new_ind,]

for(i in 1:nrow(metadata) ){
  #i=1
  print(i)
  pcm=as.matrix( read.table(paste0( metadata$filename[i]), header=FALSE) )
  
  dimnames(pcm)=list(c("A", "C", "G", "T"))
  
  min_slice_width = floor( ncol(pcm) / 3 )
  
  artificial<-list()
  for( j in min_slice_width:(ncol(pcm) - min_slice_width) ){
    #j=5
    # Get matrix versions where slices are concatenated in all 3 relative orientations
   
    artificial[paste0(metadata$ID[i], "_",j,c("_HT2", "_HH", "_TT")) ]=get_other_pfm_slice_pair_orientations(pcm, j)
    
  }
  
  for(motif in names(artificial)){
    print(motif)
    write.table(artificial[[motif]],row.names = FALSE, col.names=FALSE, quote=FALSE,file=paste0("../../PWMs_final_version2.2/artificial_motifs/", motif, ".pfm"), sep="\t")
    write.table(artificial[[motif]],row.names = FALSE, col.names=FALSE, quote=FALSE,file=paste0("../../PWMs_final_version2.2/artificial_motifs_space/", motif, ".pfm"), sep=" ")
    
    
    transfac=universalmotif::read_matrix(file=paste0("../../PWMs_final_version2.2/artificial_motifs/", motif, ".pfm"), sep="\t", header=FALSE)
    transfac@name=motif
    

    transfac_file=paste0("../../PWMs_final_version2.2/artificial_motifs_transfac/", motif,".pfm")
    universalmotif::write_transfac(transfac, file=transfac_file, overwrite = TRUE, append = FALSE)
    
    
  }
  
  artificial_all[[metadata$ID[i]]]=artificial
  
  #Write all motifs to file:
  

}

#saveRDS(artificial_all, file="RData/artificial_all.Rds")
saveRDS(artificial_all, file="RData/artificial_version2.2.Rds")

stop()

artificial_version2.2 <- readRDS("~/projects/TFBS/RProjects/TFBS/RData/artificial_version2.2.Rds")

#Do we have artifical motifs for all true motifs

files=strsplit(dir("~/projects/TFBS/PWMS_final_version2.2/artificial_motifs_transfac/"), "_")

max_length <- max(sapply(files, length))

# 2. Pad shorter elements with NA to match the longest element
padded_list <- lapply(files, function(x) {
  length(x) <- max_length
  x
})

# 3. Combine the list elements into a matrix/data frame
combined_data <- do.call(rbind, padded_list)

remove_ind=which(combined_data %in% c("HH.pfm", "HT2.pfm",  "TT.pfm"))
combined_data[remove_ind]=NA

#combined_data=apply(combined_data, 1, function(x) paste0(x, collapse="_") )


#The index of last non-na element in each row
last_non_na=apply(combined_data, 1, function(x) max(which(!is.na(x))) )
#convert the last non-na elements into na 
for(i in 1:nrow(combined_data)){
  combined_data[i,last_non_na[i]]=NA
  
}

combined_data=apply(combined_data, 1, function(x) paste0(x, collapse="_") )
combined_data <- sub("_NA$", "", combined_data)
combined_data <- sub("_NA$", "", combined_data)
combined_data <- sub("_NA$", "", combined_data)
combined_data <- sub("_NA$", "", combined_data)

length(unique(combined_data))# 3933 CORRECT

combined_data[which(!(unique(combined_data) %in% metadata$ID))]

data_path="/scratch/project_2006203/"

#df_motif_info <- read_tsv("../metadata/df_motif_info.tsv", col_names=TRUE)
#representatives <- read_tsv("../metadata/new_representatives.tsv", col_names=TRUE)
representatives <- read_tsv(paste0(data_path,"motif-clustering-Viestra-private/metadata/new_representatives_IC_length_Morgunova.tsv"), col_names=TRUE)

for(i in 1:nrow(representatives)){
  #i=1
  print(i)
  pcm=as.matrix( read.table(paste0(data_path,"TFBS/", representatives$filename[i]), header=FALSE) )
  length=ncol(pcm)
  
  #In how many different combinations you can have #length columns
  factorial(length) #40320
  
  #Generate all permutations
  library(gtools)
  perms <- gtools::permutations(ncol(pcm), ncol(pcm))
  perms <- gtools::permutations(ncol(pcm), ncol(pcm))
  pairings <- expand.grid(1:col(perms), 1:ncol(perms))
  
  all.perms <- lapply(1:nrow(pairings), function(x) pcm[, perms[pairings[x,1],] ])
  
  all.unique.perms <- unique(perms)
  length(all.unique.perms)
  
  perms2=Permn(1:length)
  perms=perms[-1,] #the original not wanted
  
  motif_perms=pcm[,perms[1,]]
  
  #For creating the whole dataset of the permutations, there’s the function DescTools::CombSet(), which has
  #an argument m for defining the size of the subset to be drawn and where the replacement and order
  #arguments can be set. 
  #CombN(1:length, 4, repl=FALSE, ord=FALSE) 

  #Sample two or three columns at a time
  
  #How many two column parts?
  
  #length(which((representatives$length %% 2) == 0)) #even 4218
  #length(which((representatives$length %% 2) != 0)) #odd 1565
  
  if(length %% 2==0){ #even
    
    parts=length/2 #length is even
    
    ind=as.matrix(data.frame(first=seq(1,(length-1),2), second=seq(2,length,2)))
    colnames(ind)=NULL
    ind.list <- split(ind, seq(nrow(ind)))
    
    ind2=as.matrix(data.frame(first=seq(2,length,2), second=c( seq(3,length-1,2),1) ))
    colnames(ind2)=NULL
    ind2.list <- split(ind2, seq(nrow(ind2)))
    
  }else{ #odd
    
    parts=ceiling(length/2)
    
    #first position alone
    ind=as.matrix(data.frame(first=c(1, seq(2,(length-1),2)), second=c(NA,seq(3,length,2))) )
    colnames(ind)=NULL
    ind.list <- split(ind, seq(nrow(ind)))
    ind.list[[1]]=1
    
    #last position alone
    ind2=as.matrix(data.frame(first=c(seq(1,length-2,2),length), second=c( seq(2,length-1,2),length) ))
    colnames(ind2)=NULL
    ind2.list <- split(ind2, seq(nrow(ind2)))
    ind2.list[[parts]]=length
    
  }
  
  #How many combinations of parts
  perms=Permn(1:parts)
  perms=perms[-1,] #the original not wanted
  
  #all possible motif permutations
  
  

  
  #reverse complements
  #A control set of reversed but not complemented motifs was generated by reversing the column order of each motif matches (Sahu et al. 2022)
  
  #handle cap-selex motifs
  

}

# x <- letters[1:4]
# n <- length(x)
# m <- 2
# factorial(n) #24
# Permn(x) 
# CombN(n, m, repl=FALSE, ord=TRUE) #12
# CombSet(x, m, repl=FALSE, ord=TRUE) #All combinations of two letters, the letter can not occur twice, the order matters
# 
# CombN(n, m, repl=TRUE, ord=TRUE) #16
# CombSet(x, m, repl=TRUE, ord=TRUE) #All combinations of two letters, the letter can not occur twice, the order matters
# 
# CombN(n, m, repl=TRUE, ord=FALSE) #10
# CombSet(x, m, repl=TRUE, ord=FALSE) #All combinations of two letters, the letter can occur twice, the order does not matter
# 
# CombN(n, m, repl=FALSE, ord=FALSE) #6
# CombSet(x, m, repl=FALSE, ord=FALSE) #All combinations of two letters, the letter can not occur twice, the order does not matter






