# =========================================================================
# 00_make_design.R
#
# Space-filling design over the community-level parameters, one row per
# parameter set. Species-level traits (K, r and mu of each species) are
# drawn at random inside each set by 02_run_batch.R, so every set keeps
# its within-community variation. This script controls how the sets are
# spread over the community-level axes:
#
#   alpha_cv    spread of the alpha matrix (SD = cv * mean)        [0.10, 0.5]
#   mean_alpha  mean competition between species                   [0.05, 1.0]
#   k_decades   width of the log10 K range, centred on the K range
#               set in the submit script (1.0 = the full range;
#               smaller values make species more similar in K)     [0.3, 1.0]
#
# To drop a dimension, set its minimum equal to its maximum (for example
# K_DECADES_MIN = K_DECADES_MAX = 1 uses the full K range in every set).
#
# Method: scrambled Sobol sequence (random digital shift) implemented in
# base R, with Joe-Kuo direction numbers for up to 5 dimensions. Any block
# of 2^m consecutive rows starting at a multiple of 2^m (1-512, 513-1024,
# ...) is itself evenly spread, so a run can be split into array batches
# along those boundaries. Coverage checks are printed at the end.
#
# Usage: Rscript 00_make_design.R
#   Environment variables override the defaults:
#   N_SETS (1024), DESIGN_SEED (2026), DESIGN_OUT (design.csv),
#   ALPHA_CV_MIN/MAX, ALPHA_MEAN_MIN/MAX, K_DECADES_MIN/MAX
# =========================================================================

N_SETS      <- as.integer(Sys.getenv("N_SETS", "1024"))
DESIGN_SEED <- as.integer(Sys.getenv("DESIGN_SEED", "2026"))
DESIGN_OUT  <- Sys.getenv("DESIGN_OUT", "design.csv")

ranges <- list(
  alpha_cv   = c(as.numeric(Sys.getenv("ALPHA_CV_MIN", "0.1")),    as.numeric(Sys.getenv("ALPHA_CV_MAX", "0.5"))),
  mean_alpha = c(as.numeric(Sys.getenv("ALPHA_MEAN_MIN", "0.05")), as.numeric(Sys.getenv("ALPHA_MEAN_MAX", "1.0"))),
  k_decades  = c(as.numeric(Sys.getenv("K_DECADES_MIN", "0.3")),   as.numeric(Sys.getenv("K_DECADES_MAX", "1.0")))
)
d <- length(ranges)

# ---- scrambled Sobol sequence in base R ---------------------------------
sobol_points <- function(n, d, seed = 1L, bits = 30L) {
  JK <- list(NULL,
             list(s = 1L, a = 0L, m = 1L),
             list(s = 2L, a = 1L, m = c(1L, 3L)),
             list(s = 3L, a = 1L, m = c(1L, 3L, 1L)),
             list(s = 3L, a = 2L, m = c(1L, 1L, 1L)))
  if (d > 5) stop("this Sobol generator supports up to 5 dimensions")
  V <- matrix(0L, nrow = bits, ncol = d)
  for (j in seq_len(d)) {
    m <- numeric(bits)
    if (j == 1) {
      m[] <- 1
    } else {
      s <- JK[[j]]$s; a <- JK[[j]]$a
      m[seq_len(s)] <- JK[[j]]$m
      for (i in (s + 1):bits) {
        mi <- bitwXor(as.integer(2^s * m[i - s]), as.integer(m[i - s]))
        for (k in seq_len(s - 1L)) {
          if (bitwAnd(bitwShiftR(a, s - 1L - k), 1L) == 1L) mi <- bitwXor(mi, as.integer(2^k * m[i - k]))
        }
        m[i] <- mi
      }
    }
    V[, j] <- as.integer(m * 2^(bits - seq_len(bits)))
  }
  pts <- matrix(0L, nrow = n, ncol = d)
  x <- integer(d)
  if (n > 1) for (k in 1:(n - 1)) {
    cpos <- 1L; kk <- k - 1L
    while (bitwAnd(kk, 1L) == 1L) { kk <- bitwShiftR(kk, 1L); cpos <- cpos + 1L }
    x <- bitwXor(x, V[cpos, ])
    pts[k + 1, ] <- x
  }
  set.seed(seed)                                   # random digital shift
  shift <- as.integer(floor(runif(d) * 2^bits))
  for (j in seq_len(d)) pts[, j] <- bitwXor(pts[, j], shift[j])
  pts / 2^bits
}

# Fall back to a Latin hypercube if there are more dimensions than the
# Sobol generator supports
if (d <= 5) {
  method <- "Sobol (scrambled, base R)"
  u <- sobol_points(N_SETS, d, seed = DESIGN_SEED)
} else {
  method <- "Latin hypercube"
  set.seed(DESIGN_SEED)
  u <- sapply(seq_len(d), function(j) (sample.int(N_SETS) - runif(N_SETS)) / N_SETS)
}
u <- matrix(u, nrow = N_SETS, ncol = d)

design <- data.frame(set_id = seq_len(N_SETS))
for (j in seq_len(d)) {
  nm <- names(ranges)[j]
  design[[nm]] <- ranges[[nm]][1] + u[, j] * (ranges[[nm]][2] - ranges[[nm]][1])
}
write.csv(design, DESIGN_OUT, row.names = FALSE)

cat("Method:", method, "\n")
cat("Wrote", N_SETS, "parameter sets to", DESIGN_OUT, "\n\n")
print(summary(design[, names(ranges)]))
cat("\nCorrelations between dimensions (should be near 0):\n")
print(round(cor(design[, names(ranges)]), 3))

# Coverage check: points per cell on an 8 x 8 grid for each pair of
# dimensions (ideal: N_SETS / 64 per cell)
cat("\nPoints per cell on an 8 x 8 grid, per pair (ideal", N_SETS / 64, "each):\n")
for (pr in combn(d, 2, simplify = FALSE)) {
  tab <- table(floor(u[, pr[1]] * 8), floor(u[, pr[2]] * 8))
  cat(sprintf("  %-10s x %-10s min %d, max %d\n", names(ranges)[pr[1]], names(ranges)[pr[2]], min(tab), max(tab)))
}

# Pairwise scatterplot of the design (skipped quietly if no graphics device)
try({
  png(sub("\\.csv$", "_check.png", DESIGN_OUT), width = 1200, height = 1200, res = 150, type = "cairo")
  pairs(design[, names(ranges)], pch = 16, cex = 0.4, gap = 0.3,
        main = paste(method, "-", N_SETS, "sets"))
  dev.off()
}, silent = TRUE)
