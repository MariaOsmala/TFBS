#Hypergeometric test/Fisher's exact test

#m: a list of all cCREs that have a motif match

#n: A set of cCREs restrictive to a given cell type. G_0, n=|G_0|

#N: G: all cCREs across all cell types G, N=|G|

#k: number of cell-line-specific cCREs that have motif match



FishersExact <- function(k,n,m,N){
  e = log2( (k / n) / (m / N) )
  
  Nm=N-m
  
  #k=oFg_hs - 1, m=oBg_hs, N-m=nBg_hs - oBg_hs n= nFg_hs
  
  #oFg = fromIntegral $ length $ filter (`S.member` motifs) fg :: Double
  #oBg = fromIntegral $ S.size motifs :: Double
  #nFg = fromIntegral $ length fg :: Double
  #nBg = fromIntegral nBg' :: Double
  
  #phyper(q/x, m, n, k)
  #phyper(k, n, N-m, m )
  
  if(e >= 0){
    #log.p: if TRUE, probabilities p are given as log(p), natural logarithm.
    
    #The logarithm log_b(x) log_10(x) can be computed from he logarithms of x and b with respect to an arbitrary base k:
    #log_b(x)=log_e(x)/log_e(b) #log_10(p)=log_e(p)/log_e(10), phyper returns log_e(p)
    p <- phyper(k - 1, n, Nm, m, lower.tail=F, log.p=T) / log(10)
    
  }else{
    # lower.tail=TRUE (default), probabilities are P[X<=k], otherwise, P[X > k].
    #k in 0, ..., n
    # Why the lower tail TRUE? Should be false?
    p <- phyper(k, n, Nm, m, lower.tail=T, log.p=T) / log(10)
  }
  
  fe<-list()
  fe$p=p
  fe$e=e
  return(fe)
  
}