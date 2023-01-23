library(rtracklayer)

#Database stuff
setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject/")

e_matrix <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/e_matrix_human_cell_type_restricted.Rds")
p_matrix <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/p_matrix_human_cell_type_restricted.Rds")

#proportion of positive and negative enrichments
table(e_matrix<=0)

#k: number of cCREs that have motif match
#n=|G_0| A set of cCREs restrictive to a given cell type. G_0,

#m=|S| A list of 300 000 regions from motif matching. S, 

#N=|G| G: all cCREs across all cell types G,  435142


e = log2( (k / n) / (m / N) )


#k=oFg_hs - 1, m=oBg_hs, N-m=nBg_hs - oBg_hs n= nFg_hs

#oFg = fromIntegral $ length $ filter (`S.member` motifs) fg :: Double
#oBg = fromIntegral $ S.size motifs :: Double
#nFg = fromIntegral $ length fg :: Double
#nBg = fromIntegral nBg' :: Double

if(e >= 0){
  p <- phyper(k - 1, m, Nm, n, lower.tail=F, log.p=T) / log(10)
  
}else{
  # lower.tail=TRUE (default), probabilities are P[X<=k], otherwise, P[X > k].
  #k in 0, ..., n
  # Why the lower tail TRUE? Should be false?
  p <- phyper(k, m, Nm, n, lower.tail=T, log.p=T) / log(10)
}
