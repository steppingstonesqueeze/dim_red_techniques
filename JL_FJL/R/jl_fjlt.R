
# R/jl_fjlt.R

suppressWarnings(suppressMessages({
  library(Matrix)
}))

# --- JL transforms ---

embed_jl_gaussian <- function(X, m, seed=123) {
  set.seed(seed)
  n <- nrow(X); D <- ncol(X)
  W <- matrix(rnorm(D*m, mean=0, sd=1.0/sqrt(m)), nrow=D, ncol=m)
  X %*% W
}

embed_jl_rademacher <- function(X, m, seed=124) {
  set.seed(seed)
  n <- nrow(X); D <- ncol(X)
  R <- matrix(sample(c(-1,1), D*m, replace=TRUE), nrow=D, ncol=m)
  R <- R / sqrt(m)
  X %*% R
}

# Achlioptas sparse JL: each entry is
#   +1/sqrt(s) with prob 1/(2s), 0 with prob 1 - 1/s, -1/sqrt(s) with prob 1/(2s)
embed_jl_achlioptas <- function(X, m, s=3, seed=125) {
  set.seed(seed)
  n <- nrow(X); D <- ncol(X)
  probs <- c(1/(2*s), 1 - 1/s, 1/(2*s))
  vals  <- c( 1/sqrt(s), 0, -1/sqrt(s))
  # generate a sparse matrix column-wise
  ij <- list()
  x  <- numeric(0)
  for (j in 1:m) {
    # sample nonzeros
    draw <- sample(vals, D, replace=TRUE, prob=probs)
    nz   <- which(draw != 0)
    ij[[j]] <- cbind(nz, rep.int(j, length(nz)))
    x <- c(x, draw[nz])
  }
  if (length(x) == 0) {
    W <- Matrix(0, nrow=D, ncol=m, sparse=TRUE)
  } else {
    IJ <- do.call(rbind, ij)
    W <- sparseMatrix(i=IJ[,1], j=IJ[,2], x=x, dims=c(D, m))
  }
  as.matrix(X %*% W)
}

# --- FJLT (Hadamard-based) ---
# Y = sqrt(D/m) * X * Dsign * H * S
# - Dsign: random ±1 diagonal
# - H: (padded) Walsh–Hadamard transform (size P, next pow2 >= D)
# - S: column sampler (m columns) from the first D columns; if padded, we sample from first D positions.

next_pow2 <- function(x) { 2^(ceiling(log2(x))) }

embed_fjlt <- function(X, m, seed=200, sample_cols="uniform") {
  # sample_cols currently only "uniform"
  set.seed(seed)
  n <- nrow(X); D <- ncol(X)
  P <- next_pow2(D)
  # step 1: random sign flip on columns
  sgn <- sample(c(-1, 1), D, replace=TRUE)
  Xs <- X
  Xs <- sweep(Xs, 2L, sgn, `*`)

  # pad to P if needed
  if (P > D) {
    Xp <- cbind(Xs, matrix(0, nrow=n, ncol=P - D))
  } else {
    Xp <- Xs
  }

  # step 2: FWHT on rows
  Xh <- fwht_matrix(Xp)  # normalized by 1/sqrt(P)

  # step 3: sample m columns (from the first D columns to avoid padded area)
  if (m > D) stop("m must be <= D")
  cols <- sample.int(D, size=m, replace=FALSE)

  # scale to preserve norms in expectation
  Y <- Xh[, cols, drop=FALSE] * sqrt(P / m)

  Y
}
