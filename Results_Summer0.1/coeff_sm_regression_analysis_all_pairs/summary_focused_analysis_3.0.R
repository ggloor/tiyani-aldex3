# =============================================================================
# Focused Summary Analysis v3.0: coefficient.sm vs CLR(std) vs TSS(std)
#
# Replaces the old custom clr_const/tss_const approach with the proper
# coefficient.sm function from ALDEx3 v1.2.0:
#   c.mu  = c(0, k)              — intercept=0, slope=k
#   c.cor = diag(c(0, gamma^2))  — no intercept uncertainty, gamma on slope
#
# Sweeps constants (-1 to 1 by 0.1) for EACH gamma value.
# Runs coefficient.sm, CLR(std), and TSS(std) for each constant.
#
# 5 unified categories: all3, coeff, clr, tss, none
# Robustness = all_three / (all_three + coeff_only + clr_only + tss_only)
#   1 = all normalizations agree, 0 = max disagreement
# Zero p-value handling: min(nonzero)/10 before -log10 (Dr. Gloor convention)
#
# Tracks:
#   1. P-value divergence (cor, MAE, RMSE) between all method pairs
#   2. Slope & intercept of "all3" (red) dots + direction split
#   3. Overall median -log10(padj) for all three methods
#
# Output folder structure:
#   pval_analysis_3.0/
#     gamma_0.1/
#       data/    — CSV
#       plots/   — per-pairwise PDFs
#     gamma_0.3/
#       ...
# =============================================================================

library(ALDEx3)

# ---- Config -----------------------------------------------------------------
datadir <- "~/Desktop/3383/0_git/tiyani-aldex3/data/"
nsample <- 32

gammas    <- c(0.1, 0.3, 0.5)
constants <- seq(-1, 1, by = 0.1)

# ---- Main output folder -----------------------------------------------------
maindir <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/pval_analysis_3.0/"

# ---- Load datasets ----------------------------------------------------------
load(paste0(datadir, "cent.Rda"))
load(paste0(datadir, "eld.Rda"))
load(paste0(datadir, "kin.Rda"))
load(paste0(datadir, "mage.Rda"))

pairs_list <- list(
  c("cent", "mage"),
  c("eld",  "kin")
)

datasets <- list(cent = cent, eld = eld, kin = kin, mage = mage)

# ---- Helper: run analysis for one constant + one comparison -----------------
run_analysis <- function(sum_coeff, sum_clr, sum_tss, k, comp) {

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

  # Significance (using original padj)
  sig_coeff <- sum_coeff$p.val.adj < 0.05
  sig_clr   <- sum_clr$p.val.adj < 0.05
  sig_tss   <- sum_tss$p.val.adj < 0.05

  # 5 unified categories
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
  n_none  <- sum(cats == "none")

  # Robustness
  n_any_sig  <- n_all3 + n_coeff + n_clr + n_tss
  robustness <- if (n_any_sig > 0) n_all3 / n_any_sig else NA

  # Summary stats row
  stats_row <- data.frame(
    constant          = k,
    comparison        = comp,
    robustness        = robustness,
    n_all3            = n_all3,
    n_coeff_only      = n_coeff,
    n_clr_only        = n_clr,
    n_tss_only        = n_tss,
    n_none            = n_none,
    median_log_coeff  = median(log_coeff),
    median_log_clr    = median(log_clr),
    median_log_tss    = median(log_tss),
    median_es_coeff   = median(abs(sum_coeff$estimate)),
    median_es_clr     = median(abs(sum_clr$estimate)),
    median_es_tss     = median(abs(sum_tss$estimate)),
    cor_coeff_clr     = cor(log_coeff, log_clr),
    cor_coeff_tss     = cor(log_coeff, log_tss),
    cor_clr_tss       = cor(log_clr, log_tss),
    mae_coeff_clr     = mean(abs(log_coeff - log_clr)),
    mae_coeff_tss     = mean(abs(log_coeff - log_tss)),
    mae_clr_tss       = mean(abs(log_clr - log_tss))
  )

  # Slope/intercept on red dots (all3) — Coeff vs CLR and Coeff vs TSS
  all3_idx       <- cats == "all3"
  red_log_coeff  <- log_coeff[all3_idx]
  red_log_clr    <- log_clr[all3_idx]
  red_log_tss    <- log_tss[all3_idx]
  red_coeff_padj <- sum_coeff$p.val.adj[all3_idx]
  red_clr_padj   <- sum_clr$p.val.adj[all3_idx]
  red_tss_padj   <- sum_tss$p.val.adj[all3_idx]

  # Direction split: Coeff vs CLR
  coeff_more_clr <- red_coeff_padj < red_clr_padj
  clr_more_coeff <- red_clr_padj < red_coeff_padj

  n_coeff_side_clr <- sum(coeff_more_clr)
  if (n_coeff_side_clr >= 2) {
    fit <- lm(red_log_clr[coeff_more_clr] ~ red_log_coeff[coeff_more_clr])
    cs_clr_slope <- unname(coef(fit)[2]); cs_clr_int <- unname(coef(fit)[1])
    cs_clr_r2 <- summary(fit)$r.squared
  } else { cs_clr_slope <- NA; cs_clr_int <- NA; cs_clr_r2 <- NA }

  n_clr_side_coeff <- sum(clr_more_coeff)
  if (n_clr_side_coeff >= 2) {
    fit <- lm(red_log_clr[clr_more_coeff] ~ red_log_coeff[clr_more_coeff])
    cl_clr_slope <- unname(coef(fit)[2]); cl_clr_int <- unname(coef(fit)[1])
    cl_clr_r2 <- summary(fit)$r.squared
  } else { cl_clr_slope <- NA; cl_clr_int <- NA; cl_clr_r2 <- NA }

  # Direction split: Coeff vs TSS
  coeff_more_tss <- red_coeff_padj < red_tss_padj
  tss_more_coeff <- red_tss_padj < red_coeff_padj

  n_coeff_side_tss <- sum(coeff_more_tss)
  if (n_coeff_side_tss >= 2) {
    fit <- lm(red_log_tss[coeff_more_tss] ~ red_log_coeff[coeff_more_tss])
    cs_tss_slope <- unname(coef(fit)[2]); cs_tss_int <- unname(coef(fit)[1])
    cs_tss_r2 <- summary(fit)$r.squared
  } else { cs_tss_slope <- NA; cs_tss_int <- NA; cs_tss_r2 <- NA }

  n_tss_side_coeff <- sum(tss_more_coeff)
  if (n_tss_side_coeff >= 2) {
    fit <- lm(red_log_tss[tss_more_coeff] ~ red_log_coeff[tss_more_coeff])
    ts_tss_slope <- unname(coef(fit)[2]); ts_tss_int <- unname(coef(fit)[1])
    ts_tss_r2 <- summary(fit)$r.squared
  } else { ts_tss_slope <- NA; ts_tss_int <- NA; ts_tss_r2 <- NA }

  slope_row <- data.frame(
    constant        = k,
    comparison      = comp,
    n_all3          = n_all3,
    # Coeff vs CLR direction
    n_coeff_side_clr  = n_coeff_side_clr,
    coeff_clr_slope   = cs_clr_slope,
    coeff_clr_int     = cs_clr_int,
    coeff_clr_r2      = cs_clr_r2,
    n_clr_side_coeff  = n_clr_side_coeff,
    clr_coeff_slope   = cl_clr_slope,
    clr_coeff_int     = cl_clr_int,
    clr_coeff_r2      = cl_clr_r2,
    # Coeff vs TSS direction
    n_coeff_side_tss  = n_coeff_side_tss,
    coeff_tss_slope   = cs_tss_slope,
    coeff_tss_int     = cs_tss_int,
    coeff_tss_r2      = cs_tss_r2,
    n_tss_side_coeff  = n_tss_side_coeff,
    tss_coeff_slope   = ts_tss_slope,
    tss_coeff_int     = ts_tss_int,
    tss_coeff_r2      = ts_tss_r2
  )

  # Red dot details
  red_names <- rownames(sum_coeff)[all3_idx]
  detail <- data.frame(
    constant   = k,
    comparison = comp,
    feature    = red_names,
    coeff_padj = sum_coeff$p.val.adj[all3_idx],
    clr_padj   = sum_clr$p.val.adj[all3_idx],
    tss_padj   = sum_tss$p.val.adj[all3_idx],
    coeff_log10 = red_log_coeff,
    clr_log10   = red_log_clr,
    tss_log10   = red_log_tss,
    coeff_est  = sum_coeff$estimate[all3_idx],
    clr_est    = sum_clr$estimate[all3_idx],
    tss_est    = sum_tss$estimate[all3_idx],
    coeff_se   = sum_coeff$std.error[all3_idx],
    clr_se     = sum_clr$std.error[all3_idx],
    tss_se     = sum_tss$std.error[all3_idx]
  )

  # Plot data
  plot_data <- list(
    log_coeff = log_coeff, log_clr = log_clr, log_tss = log_tss,
    est_coeff = sum_coeff$estimate, se_coeff = sum_coeff$std.error,
    est_clr   = sum_clr$estimate,   se_clr   = sum_clr$std.error,
    est_tss   = sum_tss$estimate,   se_tss   = sum_tss$std.error,
    cats      = cats,
    n_all3 = n_all3, n_coeff = n_coeff, n_clr = n_clr,
    n_tss = n_tss, n_none = n_none,
    sig_coeff = sum(sig_coeff), sig_clr = sum(sig_clr), sig_tss = sum(sig_tss)
  )

  list(stats = stats_row, slope = slope_row, detail = detail, plot_data = plot_data)
}

# ---- Helper: effect panel with per-panel legend -----------------------------
plot_effect_panel <- function(est, se, cats, title_text, method_sig,
                              cnt_all3, cnt_coeff, cnt_clr, cnt_tss, cnt_none) {

  plot(se[cats == "none"], est[cats == "none"],
       col = "grey75", pch = 1, cex = 0.4,
       xlim = range(se), ylim = range(est),
       xlab = "Std Error", ylab = "Effect",
       main = title_text, cex.lab = 1.2, cex.main = 0.9)

  if (any(cats == "tss"))
    points(se[cats == "tss"], est[cats == "tss"],
           col = "forestgreen", pch = 16, cex = 0.8)
  if (any(cats == "clr"))
    points(se[cats == "clr"], est[cats == "clr"],
           col = "dodgerblue", pch = 16, cex = 0.8)
  if (any(cats == "coeff"))
    points(se[cats == "coeff"], est[cats == "coeff"],
           col = "orange", pch = 16, cex = 0.8)
  if (any(cats == "all3"))
    points(se[cats == "all3"], est[cats == "all3"],
           col = "red", pch = 16, cex = 0.8)

  abline(h = 0, lwd = 1.5)
  abline(h = c(-1, 1), lty = 3, col = "grey60")

  legend("topleft",
         legend = c(paste0("All three (", cnt_all3, ")"),
                    paste0("Coeff only (", cnt_coeff, ")"),
                    paste0("CLR only (", cnt_clr, ")"),
                    paste0("TSS only (", cnt_tss, ")"),
                    paste0("Non-sig (", cnt_none, ")"),
                    paste0("Method sig: ", method_sig)),
         col = c("red", "orange", "dodgerblue", "forestgreen", "grey75", NA),
         pch = c(16, 16, 16, 16, 1, NA), cex = 0.7, bg = "white")
}

# ---- Helper: pairwise regression plots (slope, intercept, R², combined) -----
plot_pairwise <- function(slope_tbl, stats_tbl, comp, plotdir, gamma_val) {

  sub <- slope_tbl[slope_tbl$comparison == comp, ]
  sub <- sub[order(sub$constant), ]
  st  <- stats_tbl[stats_tbl$comparison == comp, ]
  st  <- st[order(st$constant), ]
  comp_label <- gsub(" ", "_", comp)

  ttl <- paste0("coeff.sm vs CLR vs TSS | ", comp, " | gamma=", gamma_val)

  pdf(paste0(plotdir, "regression_", comp_label, "_gamma", gamma_val, ".pdf"),
      width = 10, height = 7)

  # ---- Page 1: SLOPE (Coeff vs CLR) + median overlay -------------------------
  par(mar = c(5, 5, 4, 5))
  yr <- range(c(sub$coeff_clr_slope, sub$clr_coeff_slope, 1), na.rm = TRUE)
  yr <- yr + c(-0.05, 0.05) * diff(yr)

  plot(sub$constant, sub$coeff_clr_slope, type = "o", col = "orange", pch = 16, lwd = 2,
       cex = 1.2, ylim = yr, xlab = "Constant (c.mu[2])", ylab = "Slope",
       main = paste0(ttl, "\nSlope: Coeff vs CLR red-dot regression"),
       cex.lab = 1.3, cex.main = 1.0, cex.axis = 1.1)
  lines(sub$constant, sub$clr_coeff_slope, type = "o", col = "dodgerblue", pch = 17, lwd = 2, cex = 1.2)
  abline(h = 1, col = "grey30", lty = 2, lwd = 1.5)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)

  par(new = TRUE)
  med_yr <- range(c(st$median_log_coeff, st$median_log_clr, st$median_log_tss), na.rm = TRUE)
  med_yr <- med_yr + c(-0.1, 0.1) * diff(med_yr)
  plot(st$constant, st$median_log_coeff, type = "l", col = "orange",
       lty = 3, lwd = 1.8, axes = FALSE, xlab = "", ylab = "", ylim = med_yr)
  lines(st$constant, st$median_log_clr, col = "dodgerblue", lty = 3, lwd = 1.8)
  axis(4, col.axis = "grey30", cex.axis = 0.9)
  mtext(expression("Median " * -log[10](padj)), side = 4, line = 3, cex = 0.9, col = "grey30")

  legend("topleft",
         legend = c("Coeff side slope", "CLR side slope", "slope = 1",
                    "Median Coeff -log10(p)", "Median CLR -log10(p)"),
         col = c("orange", "dodgerblue", "grey30", "orange", "dodgerblue"),
         pch = c(16, 17, NA, NA, NA),
         lty = c(1, 1, 2, 3, 3), lwd = c(2, 2, 1.5, 1.8, 1.8),
         cex = 0.75, bg = "white")

  # ---- Page 2: SLOPE (Coeff vs TSS) + median overlay -------------------------
  par(mar = c(5, 5, 4, 5))
  yr <- range(c(sub$coeff_tss_slope, sub$tss_coeff_slope, 1), na.rm = TRUE)
  yr <- yr + c(-0.05, 0.05) * diff(yr)

  plot(sub$constant, sub$coeff_tss_slope, type = "o", col = "orange", pch = 16, lwd = 2,
       cex = 1.2, ylim = yr, xlab = "Constant (c.mu[2])", ylab = "Slope",
       main = paste0(ttl, "\nSlope: Coeff vs TSS red-dot regression"),
       cex.lab = 1.3, cex.main = 1.0, cex.axis = 1.1)
  lines(sub$constant, sub$tss_coeff_slope, type = "o", col = "forestgreen", pch = 17, lwd = 2, cex = 1.2)
  abline(h = 1, col = "grey30", lty = 2, lwd = 1.5)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)

  par(new = TRUE)
  plot(st$constant, st$median_log_coeff, type = "l", col = "orange",
       lty = 3, lwd = 1.8, axes = FALSE, xlab = "", ylab = "", ylim = med_yr)
  lines(st$constant, st$median_log_tss, col = "forestgreen", lty = 3, lwd = 1.8)
  axis(4, col.axis = "grey30", cex.axis = 0.9)
  mtext(expression("Median " * -log[10](padj)), side = 4, line = 3, cex = 0.9, col = "grey30")

  legend("topleft",
         legend = c("Coeff side slope", "TSS side slope", "slope = 1",
                    "Median Coeff -log10(p)", "Median TSS -log10(p)"),
         col = c("orange", "forestgreen", "grey30", "orange", "forestgreen"),
         pch = c(16, 17, NA, NA, NA),
         lty = c(1, 1, 2, 3, 3), lwd = c(2, 2, 1.5, 1.8, 1.8),
         cex = 0.75, bg = "white")

  # ---- Page 3: R² (both pairs) -----------------------------------------------
  par(mar = c(5, 5, 4, 2))
  plot(sub$constant, sub$coeff_clr_r2, type = "o", col = "dodgerblue", pch = 16, lwd = 2,
       cex = 1.2, ylim = c(0, 1.05), xlab = "Constant (c.mu[2])", ylab = expression(R^2),
       main = paste0(ttl, "\nR-squared of red-dot regression"),
       cex.lab = 1.3, cex.main = 1.0, cex.axis = 1.1)
  lines(sub$constant, sub$coeff_tss_r2, type = "o", col = "forestgreen", pch = 17, lwd = 2, cex = 1.2)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)
  legend("bottomleft",
         legend = c("Coeff vs CLR", "Coeff vs TSS"),
         col = c("dodgerblue", "forestgreen"), pch = c(16, 17),
         lty = 1, lwd = 2, cex = 0.9, bg = "white")

  # ---- Page 4: Robustness across constants ------------------------------------
  par(mar = c(5, 5, 4, 2))
  plot(st$constant, st$robustness, type = "o", col = "red", pch = 16, lwd = 2.5,
       cex = 1.2, ylim = c(0, 1.05), xlab = "Constant (c.mu[2])",
       ylab = "all_three / (all_three + coeff + clr + tss)",
       main = paste0(ttl, "\nRobustness (1 = all agree, 0 = max disagreement)"),
       cex.lab = 1.2, cex.main = 1.0, cex.axis = 1.1)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)

  # ---- Page 5: Median |effect size| across constants --------------------------
  par(mar = c(5, 5, 4, 2))
  es_yr <- range(c(st$median_es_coeff, st$median_es_clr, st$median_es_tss), na.rm = TRUE)
  es_yr <- es_yr + c(-0.05, 0.05) * diff(es_yr)

  plot(st$constant, st$median_es_coeff, type = "o", col = "orange", pch = 16, lwd = 2,
       cex = 1.2, ylim = es_yr, xlab = "Constant (c.mu[2])",
       ylab = "Median |effect size|",
       main = paste0(ttl, "\nMedian absolute effect size"),
       cex.lab = 1.3, cex.main = 1.0, cex.axis = 1.1)
  lines(st$constant, st$median_es_clr, type = "o", col = "dodgerblue", pch = 17, lwd = 2, cex = 1.2)
  lines(st$constant, st$median_es_tss, type = "o", col = "forestgreen", pch = 15, lwd = 2, cex = 1.2)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)
  legend("topleft",
         legend = c("Coeff median |est|", "CLR median |est|", "TSS median |est|"),
         col = c("orange", "dodgerblue", "forestgreen"), pch = c(16, 17, 15),
         lty = 1, lwd = 2, cex = 0.9, bg = "white")

  dev.off()
}

# ---- Helper: scatter plot PDF (one page per constant) -----------------------
plot_scatter_pdf <- function(pd_list, comp, plotdir, gamma_val) {
  comp_label <- gsub(" ", "_", comp)
  pdf(paste0(plotdir, "scatter_", comp_label, "_gamma", gamma_val, ".pdf"),
      width = 14, height = 7)

  sig_thresh <- -log10(0.05)

  for (item in pd_list) {
    if (item$comp != comp) next
    k  <- item$k
    pd <- item$pd

    cats <- pd$cats
    n <- length(cats)
    scatter_cols <- rep("grey75", n)
    scatter_cols[cats == "tss"]   <- "forestgreen"
    scatter_cols[cats == "clr"]   <- "dodgerblue"
    scatter_cols[cats == "coeff"] <- "orange"
    scatter_cols[cats == "all3"]  <- "red"
    scatter_pch <- ifelse(cats == "none", 1, 16)
    scatter_cex <- ifelse(cats == "none", 0.4, 0.7)

    robust <- if ((pd$n_all3 + pd$n_coeff + pd$n_clr + pd$n_tss) > 0) {
      round(pd$n_all3 / (pd$n_all3 + pd$n_coeff + pd$n_clr + pd$n_tss), 3)
    } else NA

    par(mfrow = c(1, 2), mar = c(5, 5, 3, 2), oma = c(0, 0, 3, 0))

    # Left: Coeff vs CLR
    xy_range <- range(c(pd$log_clr, pd$log_coeff))
    plot(pd$log_clr, pd$log_coeff,
         col = scatter_cols, pch = scatter_pch, cex = scatter_cex,
         xlim = xy_range, ylim = xy_range,
         xlab = expression(-log[10](FDR) ~ CLR),
         ylab = expression(-log[10](FDR) ~ Coefficient),
         main = "Coeff vs CLR", cex.lab = 1.2, cex.main = 0.95)
    abline(0, 1, lty = 2, lwd = 1.5)
    abline(h = sig_thresh, col = "grey60", lty = 3)
    abline(v = sig_thresh, col = "grey60", lty = 3)
    legend("topleft",
           legend = c(paste0("All three (", pd$n_all3, ")"),
                      paste0("Coeff only (", pd$n_coeff, ")"),
                      paste0("CLR only (", pd$n_clr, ")"),
                      paste0("TSS only (", pd$n_tss, ")"),
                      paste0("Non-sig (", pd$n_none, ")")),
           col = c("red", "orange", "dodgerblue", "forestgreen", "grey75"),
           pch = c(16, 16, 16, 16, 1), cex = 0.75, bg = "white")

    # Right: Coeff vs TSS
    xy_range2 <- range(c(pd$log_tss, pd$log_coeff))
    plot(pd$log_tss, pd$log_coeff,
         col = scatter_cols, pch = scatter_pch, cex = scatter_cex,
         xlim = xy_range2, ylim = xy_range2,
         xlab = expression(-log[10](FDR) ~ TSS),
         ylab = expression(-log[10](FDR) ~ Coefficient),
         main = "Coeff vs TSS", cex.lab = 1.2, cex.main = 0.95)
    abline(0, 1, lty = 2, lwd = 1.5)
    abline(h = sig_thresh, col = "grey60", lty = 3)
    abline(v = sig_thresh, col = "grey60", lty = 3)
    legend("topleft",
           legend = c(paste0("All three (", pd$n_all3, ")"),
                      paste0("Coeff only (", pd$n_coeff, ")"),
                      paste0("CLR only (", pd$n_clr, ")"),
                      paste0("TSS only (", pd$n_tss, ")"),
                      paste0("Non-sig (", pd$n_none, ")")),
           col = c("red", "orange", "dodgerblue", "forestgreen", "grey75"),
           pch = c(16, 16, 16, 16, 1), cex = 0.75, bg = "white")

    mtext(paste0(comp, "  |  k = ", k, "   gamma = ", gamma_val,
                 "   |   Robustness: ", robust),
          side = 3, outer = TRUE, cex = 1.0, font = 2)

    par(mfrow = c(1, 1))
  }
  dev.off()
}

# ---- Helper: 3-panel effect cross PDF (one page per constant) ---------------
plot_effect_pdf <- function(pd_list, comp, plotdir, gamma_val) {
  comp_label <- gsub(" ", "_", comp)
  pdf(paste0(plotdir, "effect_", comp_label, "_gamma", gamma_val, ".pdf"),
      width = 16, height = 7)

  for (item in pd_list) {
    if (item$comp != comp) next
    k  <- item$k
    pd <- item$pd
    cats <- pd$cats

    robust <- if ((pd$n_all3 + pd$n_coeff + pd$n_clr + pd$n_tss) > 0) {
      round(pd$n_all3 / (pd$n_all3 + pd$n_coeff + pd$n_clr + pd$n_tss), 3)
    } else NA

    par(mfrow = c(1, 3), mar = c(5, 5, 3, 1), oma = c(0, 0, 3, 0))

    plot_effect_panel(pd$est_coeff, pd$se_coeff, cats,
                      "Coefficient.sm effect space", pd$sig_coeff,
                      pd$n_all3, pd$n_coeff, pd$n_clr, pd$n_tss, pd$n_none)
    plot_effect_panel(pd$est_clr, pd$se_clr, cats,
                      "CLR(std) effect space", pd$sig_clr,
                      pd$n_all3, pd$n_coeff, pd$n_clr, pd$n_tss, pd$n_none)
    plot_effect_panel(pd$est_tss, pd$se_tss, cats,
                      "TSS(std) effect space", pd$sig_tss,
                      pd$n_all3, pd$n_coeff, pd$n_clr, pd$n_tss, pd$n_none)

    mtext(paste0(comp, "  |  k = ", k, "   gamma = ", gamma_val,
                 "   |   All three: ", pd$n_all3,
                 "   Robustness: ", robust),
          side = 3, outer = TRUE, cex = 1.0, font = 2)

    par(mfrow = c(1, 1))
  }
  dev.off()
}

# #############################################################################
# GAMMA SWEEP — outer loop
# #############################################################################

for (gamma_val in gammas) {

  cat("\n\n############################################################\n")
  cat("GAMMA =", gamma_val, "\n")
  cat("############################################################\n")

  # ---- Create gamma-specific subfolders ------------------------------------
  datadir_out <- paste0(maindir, "gamma_", gamma_val, "/data/")
  plotdir     <- paste0(maindir, "gamma_", gamma_val, "/plots/")
  dir.create(datadir_out, showWarnings = FALSE, recursive = TRUE)
  dir.create(plotdir,     showWarnings = FALSE, recursive = TRUE)

  # ---- Run all constants × all pairs ----------------------------------------
  summary_stats  <- data.frame()
  slope_table    <- data.frame()
  red_dot_detail <- data.frame()
  plot_data_list <- list()

  for (k in constants) {
    cat("=== constant k =", k, " ===\n")

    for (pr in pairs_list) {
      name1 <- pr[1]; name2 <- pr[2]
      comp  <- paste(name1, "vs", name2)
      cat("  ", comp, "\n")

      Y     <- cbind(datasets[[name1]], datasets[[name2]])
      conds <- c(rep(name1, ncol(datasets[[name1]])),
                 rep(name2, ncol(datasets[[name2]])))
      dat   <- data.frame(condition = conds)

      # Coefficient model
      res_coeff <- ALDEx3::aldex(Y, ~condition, dat,
                                 nsample = nsample, scale = coefficient.sm,
                                 c.mu = c(0, k),
                                 c.cor = diag(c(0, gamma_val^2)))

      # Standard CLR
      res_clr <- ALDEx3::aldex(Y, ~condition, dat,
                               nsample = nsample, scale = clr.sm, gamma = gamma_val)

      # Standard TSS
      res_tss <- ALDEx3::aldex(Y, ~condition, dat,
                               nsample = nsample, scale = tss.sm, gamma = gamma_val)

      out <- run_analysis(summary(res_coeff), summary(res_clr), summary(res_tss),
                          k, comp)
      summary_stats  <- rbind(summary_stats,  out$stats)
      slope_table    <- rbind(slope_table,    out$slope)
      red_dot_detail <- rbind(red_dot_detail, out$detail)
      plot_data_list[[length(plot_data_list) + 1]] <- list(k = k, comp = comp,
                                                            pd = out$plot_data)
    }
  }

  # ---- Save data -----------------------------------------------------------
  write.csv(summary_stats,
            paste0(datadir_out, "summary_stats_gamma", gamma_val, ".csv"),
            row.names = FALSE)
  write.csv(slope_table,
            paste0(datadir_out, "slope_table_gamma", gamma_val, ".csv"),
            row.names = FALSE)
  write.csv(red_dot_detail,
            paste0(datadir_out, "red_dot_details_gamma", gamma_val, ".csv"),
            row.names = FALSE)
  cat("Data saved to:", datadir_out, "\n")

  # ---- Generate plots per comparison ----------------------------------------
  for (comp in unique(summary_stats$comparison)) {
    plot_pairwise(slope_table, summary_stats, comp, plotdir, gamma_val)
    plot_scatter_pdf(plot_data_list, comp, plotdir, gamma_val)
    plot_effect_pdf(plot_data_list, comp, plotdir, gamma_val)
  }

  cat("Plots saved for gamma =", gamma_val, "\n")
}

cat("\n========== All done! ==========\n")
cat("Main folder:", maindir, "\n\n")
cat("Structure: gamma_{0.1,0.3,0.5}/{data,plots}/\n")
