
# R/data_generators.R

gen_gauss_iso <- function(n, D, seed=42) {
  set.seed(seed)
  matrix(rnorm(n*D), nrow=n, ncol=D)
}

gen_gauss_corr <- function(n, D, spike=10.0, seed=43) {
  # One spiked direction (axis-aligned) for simplicity
  set.seed(seed)
  X <- matrix(rnorm(n*D), nrow=n, ncol=D)
  X[,1] <- sqrt(spike) * X[,1]
  X
}

gen_sparse <- function(n, D, density=0.05, normalize=TRUE, seed=44) {
  set.seed(seed)
  X <- matrix(0, n, D)
  nnz <- max(1, floor(density * D))
  for (i in 1:n) {
    idx <- sample.int(D, nnz)
    X[i, idx] <- rnorm(nnz)
  }
  if (normalize) {
    nrms <- sqrt(rowSums(X^2))
    nrms[nrms == 0] <- 1
    X <- X / nrms
  }
  X
}
