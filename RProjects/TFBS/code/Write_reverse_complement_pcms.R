library(readr)


metadata <- read_delim("~/projects/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

metadata$revcomp_kldivergence=NA

# Function to compute KL divergence
compute_kl_divergence <- function(P, Q) {
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




for(i in 1:nrow(metadata)){
  #i=1 BCL6B
  #i=326 MAX 0.76
  #i=612 Nkx6-1 not palindromic
  #i=501 En2 is palindromic
  #which(abs(metadata$revcomp_kldivergence-3)< 0.5)
  #metadata$seed[which(metadata$revcomp_kldivergence<3)]
  #i=66 CUX1
  #i=68 CUX2
  #i=95 ERG, not palindromic
  print(i)
  pcm=as.matrix( read.table(paste0("~/projects/TFBS/", gsub("../../", "", metadata$filename[i])), header=FALSE) )
  dimnames(pcm)=list(c("A", "C", "G", "T"))
  
  # Make reverse complement
  
  # Reverse the matrix = flip the matrix horizontally
  pcm_rev=pcm[,ncol(pcm):1]
  
  #Substitute nucleotides with complements = flip the orders of A&T and C&G
  
  pcm_rev_comp=pcm_rev
  pcm_rev_comp["A", ]=pcm_rev["T",]
  pcm_rev_comp["T", ]=pcm_rev["A",]
  pcm_rev_comp["C", ]=pcm_rev["G",]
  pcm_rev_comp["G", ]=pcm_rev["C",]
    
  pcm_class_orig <- TFBSTools::PFMatrix(strand="+",
                                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25),

                                   profileMatrix=pcm
  )
  
  
  
  
  pwm_class_orig <- TFBSTools::toPWM(pcm_class_orig, type="prob", 
                                     pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))


  pcm_class_revcomp <- TFBSTools::PFMatrix(strand="+",
                                        bg=c(A=0.25, C=0.25, G=0.25, T=0.25),

                                        profileMatrix=pcm_rev_comp
  )

  
  
  
  pwm_class_revcomp <- TFBSTools::toPWM(pcm_class_revcomp, type="prob", 
                                        pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  
  
  #Compute the KL-divergence between the original and the reverse complement
  
  # Compute KL divergence
  kl_divergence1 <- compute_kl_divergence(pwm_class_orig@profileMatrix, pwm_class_revcomp@profileMatrix)
  kl_divergence2 <- compute_kl_divergence( pwm_class_revcomp@profileMatrix, pwm_class_orig@profileMatrix)
  print(paste("KL Divergence:", kl_divergence1))
  print(paste("KL Divergence:", kl_divergence2))
  
  metadata$revcomp_kldivergence[i]=kl_divergence1
  
  
  write.table(pcm_rev_comp, file=paste0("~/projects/TFBS/PWMs_final_version2.2/PWMs_revcomp_version2.2/space/",metadata$ID[i], ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
  write.table(pcm_rev_comp, file=paste0("~/projects/TFBS/PWMs_final_version2.2/PWMs_revcomp_version2.2/tab/",metadata$ID[i], ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
  
  
}

write.table(metadata, "~/projects/TFBS/PWMs_final_version2.2/metadata_representatives_revcomp_kldiv.tsv", sep="\t", quote=FALSE, row.names=FALSE)



#Maybe 3 is good threshold for palindromic motifs



library(dplyr)

hist(metadata$revcomp_kldivergence[which(metadata$new_representative=="YES")], breaks=0:263, xlim=c(0,20))

metadata$symbol[which(metadata$new_representative=="YES" & metadata$revcomp_kldivergence<3 & metadata$experiment=="HT-SELEX")]
metadata$consensus[which(metadata$new_representative=="YES" & (metadata$revcomp_kldivergence<3) & metadata$experiment=="HT-SELEX")]

which(abs(metadata$revcomp_kldivergence-4)< 0.5)
