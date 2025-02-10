




library(readr)
library(tidyverse)
library(ggplot2)

arrays <- as.numeric(commandArgs(trailingOnly = TRUE)) #1-1326

#arrays=arrays+1000

setwd("/projappl/project_2006203/TFBS/RProjects/TFBS")
#In reality 3933 motifs, typo in filename
# metadata <- read_delim("~/projects/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", 
#                        delim = "\t", escape_double = FALSE, 
#                        trim_ws = TRUE)

metadata <- read_delim("/projappl/project_2006203/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)



#Metadata for artificial half sites, 7 motifs
metadata_halfsites <- read_delim("/scratch/project_2006203/TFBS/PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/metadata.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)


#/scratch/project_2006203/TFBS/PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426
#metadata_halfsites$filename=paste0("../../", metadata_halfsites$filename) 
metadata_halfsites$filename=paste0("/scratch/project_2006203/TFBS/", metadata_halfsites$filename) 

#How to combine two data.frames that contain different column names, extract only shared columns 

#Extract shared columns

shared_columns <- intersect(names(metadata), names(metadata_halfsites))

#Combine the two data.frames

metadata <- rbind(metadata[shared_columns], metadata_halfsites[, shared_columns])





##"FOXD3" (this should be)  "IRF6"    "MAFF"    "NFATC1"  "ONECUT1" "SOX11"   "TGIF2LX"

# alternatives=list()
# alternatives[["BHLHB8"]]="BHLHA15"
# alternatives[["PROX2"]]="PROX1"
# alternatives[["NKX2-4"]]="NKX2-3"
# alternatives[["DBX2"]]="DLX2"
# alternatives[["BHLHB5"]]="BHLHA15"
# alternatives[["ZBTB33"]]="ZBTB33" #?
# alternatives[["PBX2"]]="PBX2" #?
# alternatives[["FOXI2"]]="FOXI1"
# alternatives[["gntR"]]="gntR" #?
# alternatives[["NKX2-6"]]="NKX2-3"

#print(unique(missing_monomers))
#"PBX4"   "BHLHB8" "DBX2"   "FOXO3A" "gntR"   "NKX2-4" "PBX2"   "PROX2"  "ZBTB33" "BHLHB5" "NKX2-6"




#capselex=readRDS(file="/Users/osmalama/projects/TFBS/RProjects/TFBS/RData/half_site_recognition_in_capselex_motifs_relaxed_tomtom_version2.2_artificial_halfsites.RDS")
capselex=readRDS(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/half_site_recognition_in_capselex_motifs_relaxed_tomtom_version2.2_artificial_halfsites.RDS")
composites=capselex %>% filter(study=="fromYimeng" & type=="composite") #1131
table(composites$switched_order, useNA="always")
#FALSE  TRUE  <NA> 
#  552   564    15 
table(composites$first.monomer=="")
#FALSE  TRUE 
#1116    15 
table(composites$second.monomer=="")
#FALSE 
#1131

#Is this because of the missing monomers and failure by tomtom?
#Are there still homodimers aligned against the composite motif?

table(is.na(composites$first_start)) #15
table(is.na(composites$second_start)) #15
table( is.na(composites$first_start) | is.na(composites$second_start)) #15



tmp=composites[which((is.na(composites$First_Optimal_offset) | is.na(composites$Second_Optimal_offset)) & composites$first.monomer==""),]

unique(do.call(rbind, strsplit(tmp$symbol,"_"))[,1])
#"BHLHB8" "DBX2"   "FOXO3A" "gntR"   "NKX2-4" "PBX2"   "PROX2"  "ZBTB33" These are the ones for which the monomer is missing

#For how many composite motifs, the alignment to the monomers can be defined? The alignment can not be of the spacing type.
#Answer: 1023, spacing type for 93 and no type for 15 due to missing monomer
table(composites$`type computation`)

spacing=composites %>% filter(type_computation=="spacing")

#tmp=composites %>% filter(ID=="FOXK2_ELF1_TTGAAC40NTAGC_YADIIII_NNGAAAACCGAAWNN_m1_c3b0") #Difficult composite, monomers do not align

#How many of the new composite motifs are predicted as spacing motifs


#Reorder the first.monomer and second.monomer columns if switched_order==TRUE
#capselex=capselex[, c(1,2,25:50)]

compute_kl_divergence <- function(P, Q) {
  rownames(P)=NULL
  rownames(Q)=NULL
  kl_div <- 0
  # Iterate over each position
  for (i in 1:ncol(P)) {
    # Compute KL divergence for each nucleotide at this position
    for (j in 1:nrow(P)) {
      if (P[j, i] > 0 && Q[j, i] > 0) { # Avoid log(0) and division by zero
        kl_div <- kl_div + P[j, i] * log2(P[j, i] / Q[j, i])
      }else{
        print("log(0) and division by zero")
      }
    }
  }
  return(kl_div)
}


# Function to compute KL divergence
compute_kl_divergence_another_way <- function(P, Q, kmers) {
  #Need to normalize the probabilities
  #rownames(P)=NULL
  #rownames(Q)=NULL
  kl_div <- 0
  #k=kmers[1]
  # Iterate over each position
  
  P_probs=lapply( 
    lapply(strsplit( as.character(kmers), ""), function(x) match(x, rownames(P) )), 
    function(x) P[cbind(x, 1:ncol(P))])
  
  Q_probs=lapply( 
    lapply(strsplit( as.character(kmers), ""), function(x) match(x, rownames(Q) )), 
    function(x) Q[cbind(x, 1:ncol(Q))])
  
  P_probs_prod=sapply(P_probs,prod)
  
  P_probs_log=lapply(P_probs, log2)
  Q_probs_log=lapply(Q_probs, log2)
  
  sum_P_probs_log=sapply(P_probs_log, sum)
  sum_Q_probs_log=sapply(Q_probs_log, sum)
  kl_div=sum(P_probs_prod*(sum_P_probs_log-sum_Q_probs_log))
  return(kl_div)
}

#"ATF3_ONECUT1_TCATTC40NATT_YJII_NTGACGTCANNNNATCGATN_m1_c3b0"

composites = composites %>% filter(type_computation=="composite") #1023



#save_path="~/projects/TFBS/PWMs_final_version2.2/composite_cores_and_overlapping_monomer_flanks/"
#save_path="/scratch/project_2006203/TFBS/PWMs_final_version2.2/composite_cores_and_overlapping_monomer_flanks/"
save_path="/scratch/project_2006203/TFBS/PWMs_final_version2.2/composite_cores_and_overlapping_monomer_flanks_another_way/"
#Interesting

which(composites$ID=="FOXC2_ELF2_TGACGT40NAGC_YSIII_NAGAAAACCGAAWN_m2_c3b0") #367
which(composites$ID=="PAX7_MYBL2_TATGGA40NAGA_YABI_NTAACCGATTAN_m1_c3b0u") #820
which(composites$ID=="NRL_FOSB_TTCTAA40NAAT_YWIII_WWWNTGCTGAGTCAYN_m2_c3b0u") #742
which(composites$ID=="ARX_TBX21_TCAGTG40NGAT_YYI_NAGGTGTTAATTRN_m1_c3b0") #2

caution=c()
overlap_issue=c()

#Create all possible k-mers of length 1-17

#bases <- c("A", "C", "G", "T")

# Function to generate k-mers for a given k
# generate_kmers <- function(k, bases) {
#   print(k)
#   if (k == 1) {
#     return(bases)
#   } else {
#     # Use expand.grid to create all combinations of length k
#     kmers <- do.call(paste0, expand.grid(rep(list(bases), k)))
#     return(kmers)
#   }
# }

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
#all_kmers <- lapply(1:13, generate_kmers, bases)

#saveRDS(all_kmers, file="~/projects/TFBS/RProjects/TFBS/RData/all_kmers.RDS")
#saveRDS(all_kmers, file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/all_kmers.RDS")
#all_kmers=readRDS(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/all_kmers.RDS")

#test=lapply(14, generate_kmers, bases)

#library("gtools")
# kmers=apply(permutations(n=length(bases),r=2,v=bases,repeats.allowed=T),1, paste, collapse="")
# 
# apply(kmers, 1, paste, collapse="")
# unique(t(combn(rep(bases, k), m = k)))

# [1] "20 overlap longer than monomers?"
# [1] "20: issue with overlap_end or maybe with overlap_start"
# [1] "24 overlap longer than monomers?"
# [1] "24: issue with overlap_end or maybe with overlap_start"



i=arrays
if(composites$switched_order[i]){
  first.monomer=composites$second.monomer[i]
  second.monomer=composites$first.monomer[i]
  TF1_start=composites$second_start[i]
  TF1_end=composites$second_end[i]
  TF1_orientation=composites$Second_Optimal_Orientation[i]
  TF1_length=metadata %>% filter(ID==first.monomer) %>% pull(length)
  
  TF1_start_monomer=composites$second_monomer_start[i]
  TF1_end_monomer=composites$second_monomer_end[i]
  TF2_start=composites$first_start[i]
  TF2_end=composites$first_end[i]
  TF2_orientation=composites$First_Optimal_Orientation[i]
  TF2_length=metadata %>% filter(ID==second.monomer) %>% pull(length)
  TF2_start_monomer=composites$first_monomer_start[i]
  TF2_end_monomer=composites$first_monomer_end[i]
  overlap=-(composites$overlap_length[i])
  overlap_start=composites$overlap_start[i]
  overlap_end=composites$overlap_end[i]
}else{
  first.monomer=composites$first.monomer[i]
  second.monomer=composites$second.monomer[i]
  TF1_start=composites$first_start[i]
  TF1_end=composites$first_end[i]
  TF1_orientation=composites$First_Optimal_Orientation[i]
  TF1_length=metadata %>% filter(ID==first.monomer) %>% pull(length)
  TF1_start_monomer=composites$first_monomer_start[i]
  TF1_end_monomer=composites$first_monomer_end[i]
  
  TF2_start=composites$second_start[i]
  TF2_end=composites$second_end[i]
  TF2_orientation=composites$Second_Optimal_Orientation[i]
  TF2_length=metadata %>% filter(ID==second.monomer) %>% pull(length)
  TF2_start_monomer=composites$second_monomer_start[i]
  TF2_end_monomer=composites$second_monomer_end[i]
  overlap=-(composites$overlap_length[i])
  overlap_start=composites$overlap_start[i]
  overlap_end=composites$overlap_end[i]
}

#Caution: TF1_start and TF2_start are the same or TF1_end and TF2_end are the or TF2_end is smaller then TF1_end
overlap_not_between_flanks=FALSE
if(TF1_end_monomer<TF1_length || TF2_start_monomer>1){
  caution=c(caution, i)
  overlap_not_between_flanks=TRUE
}

#The overlap can not be longer than either of the monomers
if(overlap>TF1_length || overlap>TF2_length){
  print(paste0(i, " ", "overlap longer than monomers?"))
  overlap_issue=c(overlap_issue,i)
  overlap_orig=overlap
  #9-->8
  overlap=min(overlap, TF1_length, TF2_length)
  if(TF1_end_monomer-TF1_start_monomer>overlap){
    TF1_end_monomer=TF1_end_monomer-overlap_orig+overlap
  }
  else if(TF2_end_monomer-TF2_start_monomer>overlap){
    TF2_end_monomer=TF2_end_monomer-overlap_orig+overlap
  }else{
    print("other error")
  }
  
}

#Fill the table

#read the pfm of the heterodimer

composite_core=as.matrix(read.table(composites$filename[i],sep="\t", header=FALSE))
dimnames(composite_core)=list(c("A", "C", "G", "T"))
if(length(overlap_start:overlap_end)>overlap){
  print(paste0(i,": issue with overlap_end or maybe with overlap_start"))
  overlap_end=overlap_start+overlap-1
}
composite_core=composite_core[,overlap_start:overlap_end]

if(!is.matrix(composite_core)){
  composite_core=as.matrix(composite_core, nrow=4,ncol=1)
}
#write.table(composite_core, file=paste0(save_path, composites$ID[i], "_core.pfm"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#read the pfm of the first monomer and extract and write TF1 flank

TF1_flank=as.matrix(read.table(metadata %>% filter(ID==first.monomer) %>% pull(filename),sep="\t", header=FALSE))
dimnames(TF1_flank)=list(c("A", "C", "G", "T"))
if(TF1_orientation=="+"){
  
  #Can the TF1 alignment exceed the heterodimer boundary on the right side? YES, this needs to be taken into account
  #from=TF1_length-overlap+1
  #to=TF1_length
  
  
  
}else{
  #TF1_orientation=="-"
  #Create reverse complement
  
  
  # Make reverse complement
  
  # Reverse the matrix = flip the matrix horizontally
  TF1_flank_rev=TF1_flank[,ncol(TF1_flank):1]
  
  #Substitute nucleotides with complements = flip the orders of A&T and C&G
  
  TF1_flank_rev_comp=TF1_flank_rev
  TF1_flank_rev_comp["A", ]=TF1_flank_rev["T",]
  TF1_flank_rev_comp["T", ]=TF1_flank_rev["A",]
  TF1_flank_rev_comp["C", ]=TF1_flank_rev["G",]
  TF1_flank_rev_comp["G", ]=TF1_flank_rev["C",]
  
  TF1_flank=TF1_flank_rev_comp
}

#from=TF1_length-overlap+1
#to=TF1_length

from=TF1_end_monomer-overlap+1
to=TF1_end_monomer

#The length of from:to is one, how to ensure that the result is a column vector
TF1_flank=TF1_flank[,from:to]
if(!is.matrix(TF1_flank)){
  TF1_flank=as.matrix(TF1_flank, nrow=4,ncol=1)
}

#write.table(TF1_flank, file=paste0(save_path, composites$ID[i], "_TF1_flank.pfm"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#Compute the KL-divergence between composite core and TF1 flank
pcm_class_core <- TFBSTools::PFMatrix(strand="+",
                                      bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                      profileMatrix=composite_core)
pwm_class_core <- TFBSTools::toPWM(pcm_class_core, 
                                   type="prob", 
                                   pseudocounts = 0.01, 
                                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

pcm_class_TF1_flank <- TFBSTools::PFMatrix(strand="+",
                                       bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                       profileMatrix=TF1_flank)

pwm_class_TF1_flank <- TFBSTools::toPWM(pcm_class_TF1_flank, 
                                    type="prob", 
                                    pseudocounts = 0.01, 
                                    bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

#print(ggplot() + ggseqlogo::geom_logo(  pwm_class_core@profileMatrix, method="probability", font="roboto_slab_regular", col_scheme="auto"  ) + 
#         ggseqlogo::theme_logo())
#print(ggplot() + ggseqlogo::geom_logo(  pwm_class_TF1_flank@profileMatrix, method="probability", font="roboto_slab_regular", col_scheme="auto"  ) + 
#       ggseqlogo::theme_logo())



#kl_divergence_core_flank_TF1 <- compute_kl_divergence(pwm_class_core@profileMatrix, pwm_class_TF1_flank@profileMatrix)
#kl_divergence_flank_TF1_core <- compute_kl_divergence( pwm_class_TF1_flank@profileMatrix, pwm_class_core@profileMatrix)
#kl_divergence_mean_TF1=(kl_divergence_core_flank_TF1+kl_divergence_flank_TF1_core)/2


#What is the kl-divergence against itself, the maximal possible kl-divergence?
#It is always 0, the minimal distance

#read the pfm of the second monomer and extract and write TF2 flank

TF2_flank=as.matrix(read.table(metadata %>% filter(ID==second.monomer) %>% pull(filename),sep="\t", header=FALSE))
dimnames(TF2_flank)=list(c("A", "C", "G", "T"))
if(TF2_orientation=="+"){
 
  
}else{
  
  
  # Make reverse complement
  
  # Reverse the matrix = flip the matrix horizontally
  TF2_flank_rev=TF2_flank[,ncol(TF2_flank):1]
  
  #Substitute nucleotides with complements = flip the orders of A&T and C&G
  
  TF2_flank_rev_comp=TF2_flank_rev
  TF2_flank_rev_comp["A", ]=TF2_flank_rev["T",]
  TF2_flank_rev_comp["T", ]=TF2_flank_rev["A",]
  TF2_flank_rev_comp["C", ]=TF2_flank_rev["G",]
  TF2_flank_rev_comp["G", ]=TF2_flank_rev["C",]
  
  #pcm_class_orig <- TFBSTools::PFMatrix(strand="+",
  #                                      bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
  #                                      profileMatrix=TF2_flank)
  #pwm_class_orig <- TFBSTools::toPWM(pcm_class_orig, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  
  
  #pcm_class_revcomp <- TFBSTools::PFMatrix(strand="+",
  #                                         bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
  #                                           profileMatrix=TF2_flank_rev_comp)
  #pwm_class_revcomp <- TFBSTools::toPWM(pcm_class_revcomp, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  
  #print(ggplot() + ggseqlogo::geom_logo(  pwm_class_orig@profileMatrix, method="probability", font="roboto_slab_regular", col_scheme="auto"  ) + 
  #         ggseqlogo::theme_logo())
  #print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="probability", font="roboto_slab_regular", col_scheme="auto"  ) + 
  #        ggseqlogo::theme_logo())
  
  TF2_flank=TF2_flank_rev_comp
  
}
from=1
to=overlap
TF2_flank=TF2_flank[,from:to]
if(!is.matrix(TF2_flank)){
  TF2_flank=as.matrix(TF2_flank, nrow=4,ncol=1)
}
#write.table(TF2_flank, file=paste0(save_path, composites$ID[i], "_TF2_flank.pfm"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

pcm_class_TF2_flank <- TFBSTools::PFMatrix(strand="+",
                                           bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                           profileMatrix=TF2_flank)

pwm_class_TF2_flank <- TFBSTools::toPWM(pcm_class_TF2_flank, 
                                        type="prob", 
                                        pseudocounts = 0.01, 
                                        bg=c(A=0.25, C=0.25, G=0.25, T=0.25))


#print(ggplot() + ggseqlogo::geom_logo(  pwm_class_TF2_flank@profileMatrix, method="probability", font="roboto_slab_regular", col_scheme="auto"  ) + 
#        ggseqlogo::theme_logo())








if(overlap>15){
  print(paste0(i, " too large overlap"))
  if((overlap-15)==1){
    pcm_class_core@profileMatrix=pcm_class_core@profileMatrix[,1:15]
    pcm_class_TF1_flank@profileMatrix=pcm_class_TF1_flank@profileMatrix[,1:15]
    pcm_class_TF2_flank@profileMatrix=pcm_class_TF2_flank@profileMatrix[,1:15]
    
  }
  if( ( (overlap-15)==2) ){
    pcm_class_core@profileMatrix=pcm_class_core@profileMatrix[,2:16]
    pcm_class_TF1_flank@profileMatrix=pcm_class_TF1_flank@profileMatrix[,2:16]
    pcm_class_TF2_flank@profileMatrix=pcm_class_TF2_flank@profileMatrix[,2:16]
    
  }
}

pssm_core <- TFBSTools::toPWM(pcm_class_core, type = "log2probratio", pseudocounts = 1)
pssm_TF1_flank <- TFBSTools::toPWM(pcm_class_TF1_flank, type = "log2probratio", pseudocounts = 1)
pssm_TF2_flank <- TFBSTools::toPWM(pcm_class_TF2_flank, type = "log2probratio", pseudocounts = 1)

pwm_class_TF2_flank <- TFBSTools::toPWM(pcm_class_TF2_flank, 
                                        type="prob", 
                                        pseudocounts = 0.01, 
                                        bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

pwm_class_TF1_flank <- TFBSTools::toPWM(pcm_class_TF1_flank, 
                                        type="prob", 
                                        pseudocounts = 0.01, 
                                        bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

pwm_class_core <- TFBSTools::toPWM(pcm_class_core, 
                                   type="prob", 
                                   pseudocounts = 0.01, 
                                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25))


all_kmers=readRDS(file=paste0("/scratch/project_2006203/TFBS/RProjects/TFBS/RData/",ncol(pwm_class_core@profileMatrix),"_mers.RDS"))

all_kmers=all_kmers[[1]]

#KL-divergence

kl_divergence_core_flank_TF1=compute_kl_divergence_another_way(P=pwm_class_core@profileMatrix, 
                                                               Q=pwm_class_TF1_flank@profileMatrix, 
                                                               kmers=all_kmers) 
kl_divergence_flank_TF1_core=compute_kl_divergence_another_way(P=pwm_class_TF1_flank@profileMatrix,
                                                               Q=pwm_class_core@profileMatrix, 
                                                               kmers=all_kmers) 

#the same value
#compute_kl_divergence(pwm_class_core@profileMatrix, 
#                                  pwm_class_TF1_flank@profileMatrix) 

#the same value
# compute_kl_divergence(pwm_class_TF1_flank@profileMatrix,
#                       pwm_class_core@profileMatrix) 


kl_divergence_core_flank_TF2=compute_kl_divergence_another_way(P=pwm_class_core@profileMatrix, 
                                                               Q=pwm_class_TF2_flank@profileMatrix, 
                                                               kmers=all_kmers) 
kl_divergence_flank_TF2_core=compute_kl_divergence_another_way(P=pwm_class_TF2_flank@profileMatrix,
                                                               Q=pwm_class_core@profileMatrix, 
                                                               kmers=all_kmers) 




#Shannon-Jensen divergence

mixture_core_TF1_flank=(pwm_class_core@profileMatrix+pwm_class_TF1_flank@profileMatrix)/2
mixture_core_TF2_flank=(pwm_class_core@profileMatrix+pwm_class_TF2_flank@profileMatrix)/2


kl_divergence_flank_TF1_mixture <- compute_kl_divergence_another_way(P=pwm_class_TF1_flank@profileMatrix, 
                                                                     Q=mixture_core_TF1_flank, kmers=all_kmers)
kl_divergence_core_mixture <- compute_kl_divergence_another_way(P=pwm_class_core@profileMatrix, 
                                                                Q=mixture_core_TF1_flank, kmers=all_kmers)

sh_flank_TF1_core=(kl_divergence_flank_TF1_mixture+kl_divergence_core_mixture)/2

kl_divergence_flank_TF2_mixture <- compute_kl_divergence_another_way(P=pwm_class_TF2_flank@profileMatrix, 
                                                                     Q=mixture_core_TF2_flank, kmers=all_kmers)
kl_divergence_core_mixture <- compute_kl_divergence_another_way(P=pwm_class_core@profileMatrix, 
                                                                Q=mixture_core_TF2_flank, kmers=all_kmers)

sh_flank_TF2_core=(kl_divergence_flank_TF2_mixture+kl_divergence_core_mixture)/2




kmer_scores=list()

core_kmer_scores=sapply( lapply( lapply(strsplit( as.character(all_kmers), ""), function(x) match(x, rownames(pssm_core@profileMatrix) )), 
        function(x) pssm_core@profileMatrix[cbind(x, 1:ncol(pssm_core@profileMatrix))]), sum)
names(core_kmer_scores)=all_kmers


TF1_flank_kmer_scores=sapply( lapply( lapply(strsplit( as.character(all_kmers), ""), 
                                             function(x) match(x, rownames(pssm_TF1_flank@profileMatrix) )), 
                                 function(x) pssm_TF1_flank@profileMatrix[cbind(x, 1:ncol(pssm_TF1_flank@profileMatrix))]), sum)
names(TF1_flank_kmer_scores)=all_kmers


TF2_flank_kmer_scores=sapply( lapply( lapply(strsplit( as.character(all_kmers), ""), 
                                             function(x) match(x, rownames(pssm_TF2_flank@profileMatrix) )), 
                                      function(x) pssm_TF2_flank@profileMatrix[cbind(x, 1:ncol(pssm_TF2_flank@profileMatrix))]), sum)
names(TF2_flank_kmer_scores)=all_kmers

kmer_scores[["core"]]=core_kmer_scores
kmer_scores[["TF1_flank"]]=TF1_flank_kmer_scores
kmer_scores[["TF2_flank"]]=TF2_flank_kmer_scores



composites_correct_table=data.frame(composite=composites$ID[i],
                                          first.monomer=first.monomer, 
                                          second.monomer=second.monomer, 
                                          TF1_start=TF1_start, 
                                          TF1_end=TF1_end, 
                                          TF1_orientation=TF1_orientation,
                                          TF1_length=TF1_length,
                                          TF1_start_monomer=TF1_start_monomer,
                                          TF1_end_monomer= TF1_end_monomer,
                                          TF2_start=TF2_start,
                                          TF2_end=TF2_end,
                                          TF2_orientation=TF2_orientation,
                                          TF2_length=TF2_length,
                                          TF2_start_monomer=TF2_start_monomer,
                                          TF2_end_monomer=TF2_end_monomer,
                                          overlap=overlap,
                                          overlap_start=overlap_start, 
                                          overlap_end=overlap_end,
                                          overlap_not_between_flanks=overlap_not_between_flanks,
                                          kl_divergence_core_flank_TF1=kl_divergence_core_flank_TF1,
                                          kl_divergence_flank_TF1_core=kl_divergence_flank_TF1_core,
                                          kl_divergence_core_flank_TF2=kl_divergence_core_flank_TF2,
                                          kl_divergence_flank_TF2_core=kl_divergence_flank_TF2_core,
                                          sj_divergence_flank_TF1_core=sh_flank_TF1_core,
                                          sj_divergence_flank_TF2_core=sh_flank_TF2_core,
                                          
                                          stringsAsFactors=FALSE)


  
  

save(kmer_scores,composites_correct_table, file=paste0(save_path, composites$ID[i], ".RData"))

