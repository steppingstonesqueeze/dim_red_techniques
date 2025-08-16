
# scripts/run_series_A.R
suppressWarnings(suppressMessages({
  library(optparse)
  library(ggplot2)
  library(Matrix)
  library(matrixStats)
}))

source("R/data_generators.R")
source("R/fwht.R")
source("R/jl_fjlt.R")
source("R/metrics_plots.R")

option_list <- list(
  make_option(c("--dataset"), type="character", default="iso_gauss",
              help="iso_gauss | corr_gauss | sparse"),
  make_option(c("--n"), type="integer", default=3000,
              help="number of points"),
  make_option(c("--D"), type="integer", default=512,
              help="ambient dimension"),
  make_option(c("--md"), type="character", default="0.05,0.1,0.2",
              help="comma-separated list of m/D ratios, e.g., 0.02,0.05,0.1"),
  make_option(c("--methods"), type="character", default="jl_gauss,jl_rademacher,fjlt",
              help="comma-separated: jl_gauss,jl_rademacher,jl_achli,fjlt"),
  make_option(c("--pairs"), type="integer", default=50000,
              help="number of random pairs for distance stats"),
  make_option(c("--achli_s"), type="integer", default=3,
              help="s parameter for Achlioptas sparse JL"),
  make_option(c("--seed"), type="integer", default=777,
              help="random seed master"),
  make_option(c("--outdir"), type="character", default="results",
              help="output directory for figures")
)

opt <- parse_args(OptionParser(option_list=option_list))

set.seed(opt$seed)

# --- data ---
if (opt$dataset == "iso_gauss") {
  X <- gen_gauss_iso(opt$n, opt$D, seed=opt$seed + 1)
} else if (opt$dataset == "corr_gauss") {
  X <- gen_gauss_corr(opt$n, opt$D, spike=10.0, seed=opt$seed + 2)
} else if (opt$dataset == "sparse") {
  X <- gen_sparse(opt$n, opt$D, density=0.05, normalize=TRUE, seed=opt$seed + 3)
} else {
  stop("unknown dataset")
}

# normalize rows? not required for JL; leave as-is.

md_vals <- as.numeric(strsplit(opt$md, ",")[[1]])
methods <- strsplit(opt$methods, ",")[[1]]

# Prepare distance pairs (same for all methods)
pairs <- sample_pairs(nrow(X), n_pairs=opt$pairs, seed=opt$seed + 100)
d_orig <- pairwise_dists_for_pairs(X, pairs, chunk=5000)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

for (md in md_vals) {
  m <- max(1, as.integer(md * opt$D))
  for (meth in methods) {
    cat(sprintf("\n[SeriesA] dataset=%s, D=%d, m=%d (md=%.3f), method=%s\n",
                opt$dataset, opt$D, m, md, meth))

    t0 <- proc.time()[3]
    if (meth == "jl_gauss") {
      Y <- embed_jl_gaussian(X, m, seed=opt$seed + 10 + m)
      label <- "JL (Gaussian)"
    } else if (meth == "jl_rademacher") {
      Y <- embed_jl_rademacher(X, m, seed=opt$seed + 20 + m)
      label <- "JL (Rademacher)"
    } else if (meth == "jl_achli") {
      Y <- embed_jl_achlioptas(X, m, s=opt$achli_s, seed=opt$seed + 30 + m)
      label <- sprintf("JL (Achlioptas s=%d)", opt$achli_s)
    } else if (meth == "fjlt") {
      Y <- embed_fjlt(X, m, seed=opt$seed + 40 + m)
      label <- "FJLT (Hadamard)"
    } else {
      stop("unknown method: ", meth)
    }
    t1 <- proc.time()[3]

    # distances and errors
    d_proj <- pairwise_dists_for_pairs(Y, pairs, chunk=5000)
    relerr <- relative_error(d_orig, d_proj)

    # quick stats
    stats <- c(
      md = md,
      m  = m,
      method = meth,
      median_abs_rel = median(abs(relerr)),
      p95_abs_rel = quantile(abs(relerr), 0.95),
      max_abs_rel = max(abs(relerr)),
      proj_sec = as.numeric(t1 - t0)
    )
    print(stats)

    # filenames
    tag <- sprintf("%s_%s_D%d_m%d_pairs%d", opt$dataset, meth, opt$D, m, opt$pairs)
    f1 <- file.path(opt$outdir, paste0("hist_before_after_", tag, ".png"))
    f2 <- file.path(opt$outdir, paste0("scatter_orig_vs_proj_", tag, ".png"))
    f3 <- file.path(opt$outdir, paste0("relerr_hist_", tag, ".png"))
    f4 <- file.path(opt$outdir, paste0("relerr_dot_", tag, ".png"))

    g1 <- plot_hist_before_after(d_orig, d_proj,
                                 title=sprintf("%s | %s | D=%d → m=%d (md=%.2f)",
                                               opt$dataset, label, opt$D, m, md))
    ggsave(f1, g1, width=7, height=5, dpi=140)

    g2 <- plot_scatter_orig_vs_proj(d_orig, d_proj,
                                    title=sprintf("%s | %s | D=%d → m=%d (md=%.2f)",
                                                  opt$dataset, label, opt$D, m, md))
    ggsave(f2, g2, width=7, height=5, dpi=140)

    g3 <- plot_relerr_hist(relerr,
                           title=sprintf("%s | %s | D=%d → m=%d (md=%.2f)",
                                         opt$dataset, label, opt$D, m, md))
    ggsave(f3, g3, width=7, height=5, dpi=140)

    g4 <- plot_relerr_dot(d_orig, relerr,
                          title=sprintf("%s | %s | D=%d → m=%d (md=%.2f)",
                                        opt$dataset, label, opt$D, m, md))
    ggsave(f4, g4, width=7, height=5, dpi=140)
  }
}

cat("\nDone. Check the results/ directory for PNGs.\n")
