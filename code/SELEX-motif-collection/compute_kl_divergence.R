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