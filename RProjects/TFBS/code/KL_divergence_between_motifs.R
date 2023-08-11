library(matrixStats)
library(pracma)

# KL-divergence function for two PWMs
kl_divergence <- function(pwm1, pwm2) {
  # Normalize the matrices
  pwm1 <- pwm1 / rowSums(pwm1)
  pwm2 <- pwm2 / rowSums(pwm2)
  
  # Compute the KL divergence for each position
  kl_divergences <- pwm1 * log(pwm1 / pwm2)
  
  # Compute the mean KL divergence
  average_kl_divergence <- mean(kl_divergences, na.rm = TRUE)
  
  return(average_kl_divergence)
}

# Function to create the shifted version of a matrix
shift_matrix <- function(mat, shift_rows, shift_cols) {
  num_rows <- nrow(mat)
  num_cols <- ncol(mat)
  
  # Perform the shift
  mat <- rbind(mat[(num_rows-shift_rows+1):num_rows, ], mat[1:(num_rows-shift_rows), ])
  mat <- cbind(mat[, (num_cols-shift_cols+1):num_cols], mat[, 1:(num_cols-shift_cols)])
  
  return(mat)
}

# Function to compare all possible shifts, rotations, and reflections of two PWMs
compare_transforms <- function(pwm1, pwm2) {
  num_rows <- nrow(pwm1)
  num_cols <- ncol(pwm1)
  
  best_kl <- Inf
  best_transform <- c("shift" = c(0, 0), "rotation" = 0, "reflection" = "none")
  
  # Try all possible shifts
  for (i in 0:(num_rows-1)) {
    for (j in 0:(num_cols-1)) {
      # Try all possible rotations (0, 90, 180, 270 degrees)
      for (rotation in 0:3) {
        rotated_pwm1 <- rot90(pwm1, rotation)
        # Try all possible reflections (none, horizontal, vertical)
        for (reflection in c("none", "horizontal", "vertical")) {
          if (reflection == "horizontal") {
            transformed_pwm1 <- rev(rotated_pwm1)
          } else if (reflection == "vertical") {
            transformed_pwm1 <- apply(rotated_pwm1, 2, rev)
          } else {
            transformed_pwm1 <- rotated_pwm1
          }
          # Create transformed matrix
          transformed_pwm1 <- shift_matrix(transformed_pwm1, i, j)
          
          # Compute KL divergence
          kl <- kl_divergence(transformed_pwm1, pwm2)
          
          if (kl < best_kl) {
            best_kl <- kl
            best_transform <- c("shift" = c(i, j), "rotation" = rotation, "reflection" = reflection)
          }
        }
      }
    }
  }
  
  return(list("best_kl" = best_kl, "best_transform" = best_transform))
}
