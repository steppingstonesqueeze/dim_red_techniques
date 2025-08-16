
# R/fwht.R
# In-place-like FWHT (applied row-wise) with normalization 1/sqrt(P)
# Requires number of columns to be a power of two.
fwht_matrix <- function(M) {
  n <- nrow(M); P <- ncol(M)
  if (P == 1) {
    return(M)
  }
  if (bitwAnd(P, P-1) != 0) stop("FWHT requires power-of-two number of columns")
  step <- 1L
  while (step < P) {
    jump <- step * 2L
    for (i in seq(1L, P, by=jump)) {
      idx1 <- i:(i+step-1L)
      idx2 <- (i+step):(i+jump-1L)
      A <- M[, idx1, drop=FALSE]
      B <- M[, idx2, drop=FALSE]
      M[, idx1] <- A + B
      M[, idx2] <- A - B
    }
    step <- jump
  }
  M / sqrt(P)
}
