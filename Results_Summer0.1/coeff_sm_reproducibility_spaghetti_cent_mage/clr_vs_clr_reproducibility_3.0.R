# =============================================================================
# CLR vs CLR Reproducibility Test (v3.0)
#
# Uses coefficient.sm (c.mu + c.sd) to apply constant to model coefficients.
# Compares coefficient model vs CLR(standard) vs TSS.
#
# Robustness = all_three / (all_three + coeff_only + clr_only + tss_only)
#   1 = all normalizations agree, 0 = max disagreement
#
# Plots: spaghetti lines (24 replicates) + bold control + mean ± SE ribbon
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
outdir <- paste0("/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/",
                 "Results_Summer0.1/clr_vs_clr_test_3.0/gamma", gamma_val, "/")
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

    # Coefficient model: c.mu = c(intercept, slope), c.sd = c(intercept_sd, slope_sd)
    res_coeff <- ALDEx3::aldex(Y, ~condition, dat,
                               nsample = nsample, scale = coefficient.sm,
                               c.mu = c(0, k),
                               c.sd = c(0, gamma_val))

    # Standard CLR
    res_clr <- ALDEx3::aldex(Y, ~condition, dat,
                             nsample = nsample, scale = clr.sm, gamma = gamma_val)

    # Standard TSS
    res_tss <- ALDEx3::aldex(Y, ~condition, dat,
                             nsample = nsample, scale = tss.sm, gamma = gamma_val)

    sum_coeff <- summary(res_coeff)
    sum_clr   <- summary(res_clr)
    sum_tss   <- summary(res_tss)

    # Zero p-value handling: min(nonzero)/10
    coeff_padj <- sum_coeff$p.val.adj
    clr_padj   <- sum_clr$p.val.adj
    tss_padj   <- sum_tss$p.val.adj
    all_pvals  <- c(coeff_padj, clr_padj, tss_padj)
    min_nz     <- min(all_pvals[all_pvals > 0])
    coeff_padj[coeff_padj == 0] <- min_nz / 10
    clr_padj[clr_padj == 0]     <- min_nz / 10
    tss_padj[tss_padj == 0]     <- min_nz / 10

    log_coeff <- -log10(coeff_padj)
    log_clr   <- -log10(clr_padj)
    log_tss   <- -log10(tss_padj)

    # Significance (using original padj, not zero-replaced)
    sig_coeff <- sum_coeff$p.val.adj < 0.05
    sig_clr   <- sum_clr$p.val.adj < 0.05
    sig_tss   <- sum_tss$p.val.adj < 0.05

    # 5 categories (same as effect sweep script)
    n <- length(sig_coeff)
    cats <- rep("none", n)
    cats[sig_tss & !sig_clr & !sig_coeff]   <- "tss"
    cats[sig_clr & !sig_tss & !sig_coeff]   <- "clr"
    cats[sig_coeff & !sig_clr & !sig_tss]   <- "coeff"
    cats[sig_clr & sig_tss & !sig_coeff]    <- "clr"
    cats[sig_coeff & sig_clr & !sig_tss]    <- "coeff"
    cats[sig_coeff & sig_tss & !sig_clr]    <- "coeff"
    cats[sig_coeff & sig_clr & sig_tss]     <- "all3"

    n_all3  <- sum(cats == "all3")
    n_coeff <- sum(cats == "coeff")
    n_clr   <- sum(cats == "clr")
    n_tss   <- sum(cats == "tss")

    # Robustness: 1 = all methods agree, 0 = max disagreement
    n_any_sig  <- n_all3 + n_coeff + n_clr + n_tss
    robustness <- if (n_any_sig > 0) n_all3 / n_any_sig else NA

    row <- data.frame(
      constant          = k,
      median_log_coeff  = median(log_coeff),
      median_log_clr    = median(log_clr),
      median_log_tss    = median(log_tss),
      mean_log_coeff    = mean(log_coeff),
      mean_log_clr      = mean(log_clr),
      mean_log_tss      = mean(log_tss),
      n_all3            = n_all3,
      n_coeff_only      = n_coeff,
      n_clr_only        = n_clr,
      n_tss_only        = n_tss,
      robustness        = robustness
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
csv_path <- paste0(outdir, "coeff_vs_clr_tss_replicates_gamma", gamma_val, ".csv")
write.csv(all_results, csv_path, row.names = FALSE)
cat("Data saved to:", csv_path, "\n")

# ---- Spaghetti plot helper (uses SE, not SD) --------------------------------
plot_spaghetti <- function(all_results, metric, ylab_text, title_text,
                           constants, n_replicates) {

  ctrl <- all_results[all_results$replicate == 0, ]
  ctrl <- ctrl[order(ctrl$constant), ]
  reps <- all_results[all_results$replicate > 0, ]

  # Y range
  yr <- range(all_results[[metric]], na.rm = TRUE)
  pad <- diff(yr) * 0.08
  yr <- yr + c(-pad, pad)

  # Compute ribbon (mean ± SE across replicates)
  ribbon_mat <- matrix(NA, nrow = length(constants), ncol = n_replicates)
  for (i in 1:n_replicates) {
    sub <- reps[reps$replicate == i, ]
    sub <- sub[order(sub$constant), ]
    ribbon_mat[, i] <- sub[[metric]]
  }
  ribbon_mean <- rowMeans(ribbon_mat, na.rm = TRUE)
  ribbon_sd   <- apply(ribbon_mat, 1, sd, na.rm = TRUE)
  ribbon_se   <- ribbon_sd / sqrt(n_replicates)
  ribbon_lo   <- ribbon_mean - ribbon_se
  ribbon_hi   <- ribbon_mean + ribbon_se

  # Plot
  par(mar = c(5, 5, 4, 2))
  plot(NULL, xlim = range(constants), ylim = yr,
       xlab = "Constant (c.mu[2])", ylab = ylab_text,
       main = title_text, cex.lab = 1.3, cex.main = 1.0, cex.axis = 1.1)
  grid(col = "grey85", lty = 1)

  # SE ribbon
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
                    paste0("Replicates (n=", n_replicates, ")"), "± 1 SE"),
         col = c("black", "red", rgb(0.3, 0.5, 0.8, 0.5), rgb(0.3, 0.5, 0.8, 0.15)),
         lty = c(1, 2, 1, NA), lwd = c(3, 2, 1, NA),
         pch = c(NA, NA, NA, 15), pt.cex = c(NA, NA, NA, 2),
         cex = 0.8, bg = "white")
}

# ---- Generate plots --------------------------------------------------------
pdf_path <- paste0(outdir, "coeff_vs_clr_tss_spaghetti_gamma", gamma_val, ".pdf")
pdf(pdf_path, width = 10, height = 7)

ttl_base <- paste0("coefficient.sm vs CLR(std) vs TSS - ", comp_name,
                    " (gamma=", gamma_val, ")")

# Page 1: Robustness ratio
plot_spaghetti(all_results, "robustness",
               "all_three / (all_three + coeff + clr + tss)",
               paste0(ttl_base, "\nRobustness (1 = all agree, 0 = max disagreement)"),
               constants, n_replicates)

# Page 2: Median -log10(padj) coefficient model
plot_spaghetti(all_results, "median_log_coeff",
               expression("Median " * -log[10](padj)),
               paste0(ttl_base, "\nMedian -log10(padj) coefficient model"),
               constants, n_replicates)

# Page 3: Mean -log10(padj) coefficient model
plot_spaghetti(all_results, "mean_log_coeff",
               expression("Mean " * -log[10](padj)),
               paste0(ttl_base, "\nMean -log10(padj) coefficient model"),
               constants, n_replicates)

# Page 4: n_all3 (three-way agreement count)
plot_spaghetti(all_results, "n_all3",
               "Features sig in coeff & CLR & TSS",
               paste0(ttl_base, "\nThree-way agreement count"),
               constants, n_replicates)

dev.off()
cat("Plots saved to:", pdf_path, "\n")
cat("\n========== All done! ==========\n")
