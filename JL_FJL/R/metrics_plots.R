
# R/metrics_plots.R

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(matrixStats)
}))

# Sample index pairs without replacement (approximately)
sample_pairs <- function(n, n_pairs, seed=1234) {
  set.seed(seed)
  i <- sample.int(n, n_pairs, replace=TRUE)
  j <- sample.int(n, n_pairs, replace=TRUE)
  mask <- which(i != j)
  cbind(i[mask], j[mask])
}

# Compute pairwise distances for given pairs in chunks to avoid huge memory
pairwise_dists_for_pairs <- function(X, pairs, chunk=5000) {
  n_pairs <- nrow(pairs)
  out <- numeric(n_pairs)
  start <- 1L
  while (start <= n_pairs) {
    end <- min(n_pairs, start + chunk - 1L)
    idx <- start:end
    diff <- X[pairs[idx,1], , drop=FALSE] - X[pairs[idx,2], , drop=FALSE]
    out[idx] <- sqrt(rowSums(diff * diff))
    start <- end + 1L
  }
  out
}

relative_error <- function(d_orig, d_proj) {
  # returns (proj - orig)/orig; handle zeros safely
  eps <- 1e-12
  (d_proj - d_orig) / pmax(d_orig, eps)
}

# plotting helpers ----------------------------------------------------------

plot_hist_before_after <- function(d_orig, d_proj, title="", bins=60) {
  df <- data.frame(
    dist = c(d_orig, d_proj),
    type = factor(rep(c("before", "after"), c(length(d_orig), length(d_proj))),
                  levels=c("before", "after"))
  )
  ggplot(df, aes(x=dist, fill=type)) +
    geom_histogram(alpha=0.55, position="identity", bins=bins) +
    labs(title=title, x="pairwise distance", y="count") +
    theme_minimal(base_size = 12)
}

plot_scatter_orig_vs_proj <- function(d_orig, d_proj, title="", alpha=0.2) {
  k <- length(d_orig)
  take <- if (k > 20000) sample.int(k, 20000) else 1:k
  df <- data.frame(orig=d_orig[take], proj=d_proj[take])
  ggplot(df, aes(x=orig, y=proj)) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    geom_point(alpha=alpha) +
    labs(title=title, x="original distance", y="projected distance") +
    theme_minimal(base_size = 12)
}

plot_relerr_hist <- function(relerr, title="", bins=60) {
  df <- data.frame(relerr = relerr * 100.0)
  ggplot(df, aes(x=relerr)) +
    geom_histogram(bins=bins, fill="#2b8cbe", alpha=0.8) +
    labs(title=title, x="relative error (%)", y="count") +
    theme_minimal(base_size = 12)
}

plot_relerr_dot <- function(d_orig, relerr, title="") {
  k <- length(d_orig)
  take <- if (k > 20000) sample.int(k, 20000) else 1:k
  df <- data.frame(orig=d_orig[take], rel=relerr[take] * 100.0)
  ggplot(df, aes(x=orig, y=rel)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_point(alpha=0.2) +
    labs(title=title, x="original distance", y="relative error (%)") +
    theme_minimal(base_size = 12)
}
