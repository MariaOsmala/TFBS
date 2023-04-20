library("rtracklayer")
library("dbplyr")
library("dplyr")
library("motifStack")
library("universalmotif")

#Database stuff
setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")



data_path="/scratch/project_2006203/TFBS/"

representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)

representatives$IC=c()
representatives$length=c()

pfm_list<-list()

for(i in 1:nrow(representatives)){
  pcm=as.matrix(read.csv(paste0("../../",representatives$filename[i]), header=FALSE, sep="\t"))
  colnames(pcm)=NULL
  rownames(pcm)=c("A", "C", "G", "T")
  
  #Matrix class: Zipper type?
  pfm <- TFBSTools::PFMatrix(ID=representatives$ID[i], name=representatives$symbol[i], 
                             bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                             tags=list(family=representatives$Lambert2018.families[i],
                                       study=representatives$study[i],
                                       experiment=representatives$experiment[i],
                                       organism=representatives$organism[i],
                                       representative=representatives$new_representative[i]),
                             profileMatrix=pcm
  )
  
  ## convert a PFM to PWM, #Pseudo-count adjusted log-odds scoring matrix
  #pwm <- TFBSTools::toPWM(pfm, type="log2probratio", pseudocounts=0.01)
  #[,1]       [,2]      [,3]       [,4]       [,5]       [,6]       [,7]        [,8]
  #A -1.4507888 -1.4270147  1.533854   1.937764 -19.408730 -3.8617812  1.6195158  0.05823456
  #C  0.9679349  0.5024802 -1.194585  -2.858568  -5.407234 -0.6344636 -2.1913814  0.45614799
  #G -0.1755438 -1.2526997 -1.349047  -5.011165 -19.408730 -2.0811106 -0.9908972 -0.09012369
  #T -0.3351416  0.8414026 -1.862740 -19.462440   1.991475  1.6091418 -2.2850757 -0.62706367
  
  icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
  #[,1]       [,2]       [,3]         [,4]         [,5]       [,6]       [,7]        [,8]
  #A 0.02152601 0.02703334 0.52079499 1.646003e+00 6.996658e-07 0.01609050 0.67774730 0.012736730
  #C 0.11509908 0.10297588 0.07858232 5.923671e-02 1.147522e-02 0.15069157 0.04829180 0.016781927
  #G 0.05210153 0.03050519 0.07060368 1.332277e-02 6.996658e-07 0.05528468 0.11098267 0.011492055
  #T 0.04664518 0.13024509 0.04945274 5.947401e-07 1.936396e+00 0.71364379 0.04525521 0.007920683
  
  representatives[i, "IC"]=sum(rowSums(icm)) #6.77892
  
  #motif length
  representatives[i, "length"]=length(pfm)
  pfm_list[[representatives$ID[i]]]
}

write.table(representatives, "/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length.tsv", sep="\t", row.names=FALSE)
