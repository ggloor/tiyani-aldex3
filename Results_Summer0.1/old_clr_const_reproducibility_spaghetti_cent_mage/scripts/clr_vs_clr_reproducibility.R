# =============================================================================
# CLR vs CLR Reproducibility Test
#
# Compares CLR(const=k) vs CLR(standard) across constants -1 to 1.
# Runs 1 control + 24 replicates to measure Monte Carlo variance.
#
# Metrics per constant per replicate:
#   - Median -log10(padj) for CLR(const) and CLR(standard)
#   - Mean   -log10(padj) for CLR(const) and CLR(standard)
#   - Robustness = both_sig / (both_sig + const_only)
#
# Plots: spaghetti lines (24 replicates) + bold control + mean ± SD ribbon
# =============================================================================

library(ALDEx3)

# ---- Config (easy to change) -----------------------------------------------
datadir      <- "~/Desktop/3383/0_git/tiyani-aldex3/data/"
nsample      <- 32
gamma_val    <- 0.5
constants    <- seq(-1, 1, by = 0.1)
n_replicates <- 24
comp_name    <- "cent_vs_mage"

# ---- Output ----------------------------------------------------------------
outdir <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/clr_vs_clr_test/gamma0.5/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data -------------------------------------------------------------
load(paste0(datadir, "cent.Rda"))
load(paste0(datadir, "mage.Rda"))

Y     <- cbind(cent, mage)
conds <- c(rep("cent", ncol(cent)), rep("mage", ncol(mage)))
dat   <- data.frame(condition = conds)

# ---- Helper: one full sweep across constants -------------------------------
run_one_sweep <- function(Y, dat, constants, gamma_val, nsample) {
  results <- data.frame()

  for (k in constants) {
    cat("    constant =", k, "\n")

    # CLR with constant (defined inside loop to capture k via closure)
    clr_const <- function(X, logComp, gamma = 0.5) {
      P  <- nrow(X); ns <- dim(logComp)[3]
      LambdaScale <- matrix(rnorm(P * ns, k, gamma), P, ns)
      logScale <- t(X) %*% LambdaScale
      return(logScale)
    }

    # Run both
    res_const <- ALDEx3::aldex(Y, ~condition, dat,
                               nsample = nsample, scale = clr_const, gamma = gamma_val)
    res_std   <- ALDEx3::aldex(Y, ~condition, dat,
                               nsample = nsample, scale = clr.sm, gamma = gamma_val)

    sum_const <- summary(res_const)
    sum_std   <- summary(res_std)

    # Zero p-value handling: min(nonzero)/10
    const_padj <- sum_const$p.val.adj
    std_padj   <- sum_std$p.val.adj
    all_pvals  <- c(const_padj, std_padj)
    min_nz     <- min(all_pvals[all_pvals > 0])
    const_padj[const_padj == 0] <- min_nz / 10
    std_padj[std_padj == 0]     <- min_nz / 10

    log_const <- -log10(const_padj)
    log_std   <- -log10(std_padj)

    # Significance categories (using original padj, not zero-replaced)
    sig_const  <- sum_const$p.val.adj < 0.05
    sig_std    <- sum_std$p.val.adj < 0.05
    both       <- sig_const & sig_std
    const_only <- sig_const & !sig_std
    std_only   <- !sig_const & sig_std

    n_both       <- sum(both)
    n_const_only <- sum(const_only)
    n_std_only   <- sum(std_only)
    robustness   <- if ((n_both + n_const_only) > 0) n_both / (n_both + n_const_only) else NA

    row <- data.frame(
      constant         = k,
      median_log_const = median(log_const),
      median_log_std   = median(log_std),
      mean_log_const   = mean(log_const),
      mean_log_std     = mean(log_std),
      n_both           = n_both,
      n_const_only     = n_const_only,
      n_std_only       = n_std_only,
      robustness       = robustness
    )
    results <- rbind(results, row)
  }
  return(results)
}

# ---- Run control (replicate 0) ---------------------------------------------
cat("========== Control run ==========\n")
control <- run_one_sweep(Y, dat, constants, gamma_val, nsample)
control$replicate <- 0

# ---- Run replicates --------------------------------------------------------
all_results <- control
for (r in 1:n_replicates) {
  cat("========== Replicate", r, "/", n_replicates, " ==========\n")
  rep_result <- run_one_sweep(Y, dat, constants, gamma_val, nsample)
  rep_result$replicate <- r
  all_results <- rbind(all_results, rep_result)
}

# ---- Save data -------------------------------------------------------------
csv_path <- paste0(outdir, "clr_vs_clr_replicates_gamma", gamma_val, ".csv")
write.csv(all_results, csv_path, row.names = FALSE)
cat("Data saved to:", csv_path, "\n")

# ---- Spaghetti plot helper -------------------------------------------------
plot_spaghetti <- function(all_results, metric, ylab_text, title_text,
                           constants, n_replicates) {

  ctrl <- all_results[all_results$replicate == 0, ]
  ctrl <- ctrl[order(ctrl$constant), ]
  reps <- all_results[all_results$replicate > 0, ]

  # Y range
  yr <- range(all_results[[metric]], na.rm = TRUE)
  pad <- diff(yr) * 0.08
  yr <- yr + c(-pad, pad)

  # Compute ribbon (mean ± SD across replicates)
  ribbon_mat <- matrix(NA, nrow = length(constants), ncol = n_replicates)
  for (i in 1:n_replicates) {
    sub <- reps[reps$replicate == i, ]
    sub <- sub[order(sub$constant), ]
    ribbon_mat[, i] <- sub[[metric]]
  }
  ribbon_mean <- rowMeans(ribbon_mat, na.rm = TRUE)
  ribbon_sd   <- apply(ribbon_mat, 1, sd, na.rm = TRUE)
  ribbon_lo   <- ribbon_mean - ribbon_sd
  ribbon_hi   <- ribbon_mean + ribbon_sd

  # Plot
  par(mar = c(5, 5, 4, 2))
  plot(NULL, xlim = range(constants), ylim = yr,
       xlab = "Constant", ylab = ylab_text,
       main = title_text, cex.lab = 1.3, cex.main = 1.0, cex.axis = 1.1)
  grid(col = "grey85", lty = 1)

  # SD ribbon
  polygon(c(constants, rev(constants)),
          c(ribbon_lo, rev(ribbon_hi)),
          col = rgb(0.3, 0.5, 0.8, 0.15), border = NA)

  # Replicate lines (thin, semi-transparent blue)
  for (i in 1:n_replicates) {
    sub <- reps[reps$replicate == i, ]
    sub <- sub[order(sub$constant), ]
    lines(sub$constant, sub[[metric]], col = rgb(0.3, 0.5, 0.8, 0.35), lwd = 0.8)
  }

  # Replicate mean (dashed red)
  lines(constants, ribbon_mean, col = "red", lwd = 2, lty = 2)

  # Control (bold black, on top)
  lines(ctrl$constant, ctrl[[metric]], col = "black", lwd = 3)

  abline(v = 0, col = "grey50", lty = 3)

  legend("topleft",
         legend = c("Control (run 1)", "Replicate mean",
                    paste0("Replicates (n=", n_replicates, ")"), "± 1 SD"),
         col = c("black", "red", rgb(0.3, 0.5, 0.8, 0.5), rgb(0.3, 0.5, 0.8, 0.15)),
         lty = c(1, 2, 1, NA), lwd = c(3, 2, 1, NA),
         pch = c(NA, NA, NA, 15), pt.cex = c(NA, NA, NA, 2),
         cex = 0.8, bg = "white")
}

# ---- Generate plots --------------------------------------------------------
pdf_path <- paste0(outdir, "clr_vs_clr_spaghetti_gamma", gamma_val, ".pdf")
pdf(pdf_path, width = 10, height = 7)

ttl_base <- paste0("CLR(const) vs CLR(std) - ", comp_name, " (gamma=", gamma_val, ")")

# Page 1: Robustness ratio
plot_spaghetti(all_results, "robustness",
               "both_sig / (both_sig + const_only)",
               paste0(ttl_base, "\nRobustness ratio"),
               constants, n_replicates)

# Page 2: Median -log10(padj) CLR(const)
plot_spaghetti(all_results, "median_log_const",
               expression("Median " * -log[10](padj)),
               paste0(ttl_base, "\nMedian -log10(padj) CLR(const)"),
               constants, n_replicates)

# Page 3: Mean -log10(padj) CLR(const)
plot_spaghetti(all_results, "mean_log_const",
               expression("Mean " * -log[10](padj)),
               paste0(ttl_base, "\nMean -log10(padj) CLR(const)"),
               constants, n_replicates)

# Page 4: n_both (context for robustness ratio)
plot_spaghetti(all_results, "n_both",
               "Number of features sig in both",
               paste0(ttl_base, "\nBoth-significant count"),
               constants, n_replicates)

dev.off()
cat("Plots saved to:", pdf_path, "\n")
cat("\n========== All done! ==========\n")
