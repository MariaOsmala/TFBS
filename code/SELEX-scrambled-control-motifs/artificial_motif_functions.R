reverseComplement <- function(mat){
  rc_matrix <- matrix(nrow = 4, ncol = ncol(mat) )
  
  for (i in 1:4) {
    for (j in 1:ncol(mat)) {
      rc_matrix[4 - i + 1, ncol(mat) - j + 1] <- mat[i, j]
    }
  }
  rc_matrix
}

get_other_pfm_slice_pair_orientations<- function(pcm, cut_point){
  #cut_point=j
  slice1 = pcm[,1:cut_point]; #From start to cut_point
  slice2 = pcm[,(cut_point + 1):ncol(pcm) ]; #From cut_point to end
  
  # reverse complement slices
  slice1_rc = reverseComplement(slice1)
  slice2_rc = reverseComplement(slice2)    

  # paste in all orders
  # HT1 ++, first, second, this would be the original
  #my $ht1 = paste($slice1, $slice2);
  # HT2 ++, second, first 
  ht2 = cbind(slice2, slice1)

  # HH +-, first, second rc
  hh = cbind(slice1, slice2_rc)
  # TT -+, first rc, second 
  tt = cbind(slice1_rc, slice2)
  list(ht2, hh, tt)
  
}