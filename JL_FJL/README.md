
# Series A — Johnson–Lindenstrauss (JL) & Fast JL (FJLT) in R

This repo runs JL (Gaussian, Rademacher, Achlioptas sparse) and FJLT (Hadamard-based, with auto-padding to next power-of-two) on the **same datasets** and produces:
- Before/after **histograms of pairwise distances**
- **Relative % error** dot plot and histogram
- **Faceted plots** over different target dimension ratios m/D

## Quick start

```bash
Rscript scripts/run_series_A.R --dataset iso_gauss --n 3000 --D 512 --md 0.02,0.05,0.1,0.2 --methods jl_gauss,jl_rademacher,fjlt --pairs 50000
```

Outputs (PNGs) are written to `results/`.

### Dependencies
- base R (≥ 4.0)
- packages: `ggplot2`, `Matrix`, `matrixStats`, `optparse`
  ```r
  install.packages(c("ggplot2", "Matrix", "matrixStats", "optparse"))
  ```

---

## Notes
- Pairwise distances are computed on a sampled set of pairs (default `--pairs 50k`) for scalability.
- FJLT uses a fast Walsh–Hadamard transform (FWHT). If D is not a power of two, we **pad** to the next power of two internally.
- Sparse JL uses Achlioptas' 3-sparse distribution (`s=3` default). You can change with `--achli_s`.
- Results include: histograms (before/after), scatter (d_orig vs d_proj), relative error hist, and faceted variants by m/D.
