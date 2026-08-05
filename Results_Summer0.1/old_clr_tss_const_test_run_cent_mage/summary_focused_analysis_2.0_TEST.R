# =============================================================================
# Focused Summary Analysis: cent vs mage & eld vs kin
#
# PART 1 — CLR constant: CLR(const) vs TSS(standard)
# PART 2 — TSS constant: CLR(standard) vs TSS(const)
#
# Sweeps across constants (-1 to 1 by 0.1) for EACH gamma value.
# Tracks:
#   1. P-value divergence between CLR and TSS
#   2. Slope & intercept of "Both sig" (red) dots + direction split
#   3. Overall median -log10(padj) for CLR and TSS (all features)
#
# Output folder structure:
#   focused_analysis/
#     CLR_constant/
#       gamma_0.1/
#         data/    — Rda + CSV
#         plots/   — per-pairwise PDFs
#       gamma_0.3/
#         ...
#     TSS_constant/
#       gamma_0.1/
#         ...
# =============================================================================

library(ALDEx3)

# ---- Config -----------------------------------------------------------------
datadir <- "~/Desktop/3383/0_git/tiyani-aldex3/data/"
nsample <- 32

gammas    <- c(0.3)
constants <- c(0)

# ---- Main output folder -----------------------------------------------------
maindir <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/pval_analysis0.2_TEST/"

# ---- Load datasets ----------------------------------------------------------
load(paste0(datadir, "cent.Rda"))
load(paste0(datadir, "eld.Rda"))
load(paste0(datadir, "kin.Rda"))
load(paste0(datadir, "mage.Rda"))

pairs_list <- list(
  c("cent", "mage")
)

datasets <- list(cent = cent, eld = eld, kin = kin, mage = mage)

# ---- Helper: shared analysis block ------------------------------------------
run_analysis <- function(sum_clr, sum_tss, k, comp) {

  clr_padj <- sum_clr$p.val.adj
  tss_padj <- sum_tss$p.val.adj
  all_pvals <- c(clr_padj, tss_padj)
  min_nz    <- min(all_pvals[all_pvals > 0])
  clr_padj[clr_padj == 0] <- min_nz / 10
  tss_padj[tss_padj == 0] <- min_nz / 10

  log_clr <- -log10(clr_padj)
  log_tss <- -log10(tss_padj)

  sig_clr <- sum_clr$p.val.adj < 0.05
  sig_tss <- sum_tss$p.val.adj < 0.05
  both     <- sig_clr & sig_tss
  tss_only <- sig_tss & !sig_clr
  clr_only <- sig_clr & !sig_tss
  neither  <- !sig_clr & !sig_tss

  # Summary stats row (includes overall medians)
  stats_row <- data.frame(
    constant         = k,
    comparison       = comp,
    pval_cor         = cor(log_clr, log_tss),
    pval_mae         = mean(abs(log_clr - log_tss)),
    pval_rmse        = sqrt(mean((log_clr - log_tss)^2)),
    pval_median_diff = median(log_clr - log_tss),
    median_log_clr   = median(log_clr),
    median_log_tss   = median(log_tss),
    median_es_clr    = median(abs(sum_clr$estimate)),
    median_es_tss    = median(abs(sum_tss$estimate)),
    n_both           = sum(both),
    n_tss_only       = sum(tss_only),
    n_clr_only       = sum(clr_only),
    n_neither        = sum(neither)
  )

  # Slope/intercept split on red dots only
  red_log_tss  <- log_tss[both]
  red_log_clr  <- log_clr[both]
  red_tss_padj <- sum_tss$p.val.adj[both]
  red_clr_padj <- sum_clr$p.val.adj[both]

  tss_side <- red_tss_padj < red_clr_padj
  clr_side <- red_clr_padj < red_tss_padj

  n_tss_side <- sum(tss_side)
  if (n_tss_side >= 2) {
    fit <- lm(red_log_clr[tss_side] ~ red_log_tss[tss_side])
    ts_slope <- unname(coef(fit)[2]); ts_int <- unname(coef(fit)[1]); ts_r2 <- summary(fit)$r.squared
  } else { ts_slope <- NA; ts_int <- NA; ts_r2 <- NA }

  n_clr_side <- sum(clr_side)
  if (n_clr_side >= 2) {
    fit <- lm(red_log_clr[clr_side] ~ red_log_tss[clr_side])
    cs_slope <- unname(coef(fit)[2]); cs_int <- unname(coef(fit)[1]); cs_r2 <- summary(fit)$r.squared
  } else { cs_slope <- NA; cs_int <- NA; cs_r2 <- NA }

  slope_row <- data.frame(
    constant      = k,      comparison    = comp,
    n_red_total   = sum(both),
    n_tss_side    = n_tss_side, tss_slope = ts_slope, tss_intercept = ts_int, tss_r2 = ts_r2,
    n_clr_side    = n_clr_side, clr_slope = cs_slope, clr_intercept = cs_int, clr_r2 = cs_r2,
    n_equal       = sum(red_tss_padj == red_clr_padj)
  )

  # Red dot details
  red_idx   <- which(both)
  red_names <- rownames(sum_clr)[red_idx]
  detail <- data.frame(
    constant  = k, comparison = comp, feature = red_names,
    clr_padj  = sum_clr$p.val.adj[red_idx], tss_padj  = sum_tss$p.val.adj[red_idx],
    clr_log10 = log_clr[red_idx],           tss_log10 = log_tss[red_idx],
    clr_est   = sum_clr$estimate[red_idx],  tss_est   = sum_tss$estimate[red_idx],
    clr_se    = sum_clr$std.error[red_idx], tss_se    = sum_tss$std.error[red_idx],
    direction = ifelse(red_tss_padj < red_clr_padj, "TSS_more_sig",
                ifelse(red_clr_padj < red_tss_padj, "CLR_more_sig", "equal"))
  )

  # Data for scatter and effect cross plots
  plot_data <- list(
    log_clr  = log_clr,
    log_tss  = log_tss,
    est_clr  = sum_clr$estimate,
    se_clr   = sum_clr$std.error,
    est_tss  = sum_tss$estimate,
    se_tss   = sum_tss$std.error,
    both     = both,
    tss_only = tss_only,
    clr_only = clr_only,
    neither  = neither
  )

  list(stats = stats_row, slope = slope_row, detail = detail, plot_data = plot_data)
}

# ---- Helper: save Rda split by comparison and direction ----------------------
save_split_rda <- function(slope_tbl, prefix, datadir, gamma_val) {
  cols_common <- c("constant", "n_red_total")
  cols_tss    <- c(cols_common, "n_tss_side", "tss_slope", "tss_intercept", "tss_r2")
  cols_clr    <- c(cols_common, "n_clr_side", "clr_slope", "clr_intercept", "clr_r2")
  clean_names <- c("constant", "n_red_total", "n", "slope", "intercept", "r2")

  for (comp in unique(slope_tbl$comparison)) {
    comp_label <- gsub(" ", "_", comp)

    tss_sub <- slope_tbl[slope_tbl$comparison == comp, cols_tss]
    colnames(tss_sub) <- clean_names
    rownames(tss_sub) <- NULL
    save(tss_sub, file = paste0(datadir, prefix, "_", comp_label,
                                "_TSS_side_gamma", gamma_val, ".Rda"))

    clr_sub <- slope_tbl[slope_tbl$comparison == comp, cols_clr]
    colnames(clr_sub) <- clean_names
    rownames(clr_sub) <- NULL
    save(clr_sub, file = paste0(datadir, prefix, "_", comp_label,
                                "_CLR_side_gamma", gamma_val, ".Rda"))
  }
}

# ---- Helper: plot per-pairwise PDF with median overlay -----------------------
plot_pairwise <- function(slope_tbl, stats_tbl, comp, label, plotdir, gamma_val) {

  sub  <- slope_tbl[slope_tbl$comparison == comp, ]
  sub  <- sub[order(sub$constant), ]
  st   <- stats_tbl[stats_tbl$comparison == comp, ]
  st   <- st[order(st$constant), ]
  comp_label <- gsub(" ", "_", comp)

  ttl <- paste0(label, " | ", comp, " | gamma=", gamma_val)

  pdf(paste0(plotdir, label, "_", comp_label, "_gamma", gamma_val, ".pdf"),
      width = 10, height = 7)

  # ---- Page 1: SLOPE + median overlay ----------------------------------------
  par(mar = c(5, 5, 4, 5))
  yr <- range(c(sub$tss_slope, sub$clr_slope, 1), na.rm = TRUE)
  yr <- yr + c(-0.05, 0.05) * diff(yr)

  plot(sub$constant, sub$tss_slope, type = "o", col = "blue", pch = 16, lwd = 2,
       cex = 1.2, ylim = yr, xlab = "Constant", ylab = "Slope",
       main = paste0(ttl, "\nSlope of red-dot regression"),
       cex.lab = 1.3, cex.main = 1.1, cex.axis = 1.1)
  lines(sub$constant, sub$clr_slope, type = "o", col = "red", pch = 17, lwd = 2, cex = 1.2)
  abline(h = 1, col = "grey30", lty = 2, lwd = 1.5)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)

  # Median overlay on right axis
  par(new = TRUE)
  med_yr <- range(c(st$median_log_clr, st$median_log_tss), na.rm = TRUE)
  med_yr <- med_yr + c(-0.1, 0.1) * diff(med_yr)
  plot(st$constant, st$median_log_tss, type = "l", col = "blue",
       lty = 3, lwd = 1.8, axes = FALSE, xlab = "", ylab = "", ylim = med_yr)
  lines(st$constant, st$median_log_clr, col = "red", lty = 3, lwd = 1.8)
  axis(4, col.axis = "grey30", cex.axis = 0.9)
  mtext(expression("Median " * -log[10](padj)), side = 4, line = 3, cex = 0.9, col = "grey30")

  legend("topleft",
         legend = c("TSS side slope", "CLR side slope", "slope = 1",
                    "Median TSS -log10(padj)", "Median CLR -log10(padj)"),
         col = c("blue", "red", "grey30", "blue", "red"),
         pch = c(16, 17, NA, NA, NA),
         lty = c(1, 1, 2, 3, 3), lwd = c(2, 2, 1.5, 1.8, 1.8),
         cex = 0.75, bg = "white")

  # ---- Page 2: INTERCEPT + median overlay ------------------------------------
  par(mar = c(5, 5, 4, 5))
  yr <- range(c(sub$tss_intercept, sub$clr_intercept, 0), na.rm = TRUE)
  yr <- yr + c(-0.05, 0.05) * diff(yr)

  plot(sub$constant, sub$tss_intercept, type = "o", col = "blue", pch = 16, lwd = 2,
       cex = 1.2, ylim = yr, xlab = "Constant", ylab = "Intercept",
       main = paste0(ttl, "\nIntercept of red-dot regression"),
       cex.lab = 1.3, cex.main = 1.1, cex.axis = 1.1)
  lines(sub$constant, sub$clr_intercept, type = "o", col = "red", pch = 17, lwd = 2, cex = 1.2)
  abline(h = 0, col = "grey30", lty = 2, lwd = 1.5)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)

  par(new = TRUE)
  plot(st$constant, st$median_log_tss, type = "l", col = "blue",
       lty = 3, lwd = 1.8, axes = FALSE, xlab = "", ylab = "", ylim = med_yr)
  lines(st$constant, st$median_log_clr, col = "red", lty = 3, lwd = 1.8)
  axis(4, col.axis = "grey30", cex.axis = 0.9)
  mtext(expression("Median " * -log[10](padj)), side = 4, line = 3, cex = 0.9, col = "grey30")

  legend("topleft",
         legend = c("TSS side intercept", "CLR side intercept", "intercept = 0",
                    "Median TSS -log10(padj)", "Median CLR -log10(padj)"),
         col = c("blue", "red", "grey30", "blue", "red"),
         pch = c(16, 17, NA, NA, NA),
         lty = c(1, 1, 2, 3, 3), lwd = c(2, 2, 1.5, 1.8, 1.8),
         cex = 0.75, bg = "white")

  # ---- Page 3: R² -----------------------------------------------------------
  par(mar = c(5, 5, 4, 2))

  plot(sub$constant, sub$tss_r2, type = "o", col = "blue", pch = 16, lwd = 2,
       cex = 1.2, ylim = c(0, 1.05), xlab = "Constant", ylab = expression(R^2),
       main = paste0(ttl, "\nR-squared of red-dot regression"),
       cex.lab = 1.3, cex.main = 1.1, cex.axis = 1.1)
  lines(sub$constant, sub$clr_r2, type = "o", col = "red", pch = 17, lwd = 2, cex = 1.2)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)
  legend("bottomleft",
         legend = c("TSS side", "CLR side"),
         col = c("blue", "red"), pch = c(16, 17),
         lty = 1, lwd = 2, cex = 0.9, bg = "white")

  # ---- Page 4: SLOPE + INTERCEPT combined (2-panel) + median overlay ---------
  par(mfrow = c(1, 2))

  # Left: slope + median
  par(mar = c(5, 5, 4, 5))
  yr_s <- range(c(sub$tss_slope, sub$clr_slope, 1), na.rm = TRUE)
  yr_s <- yr_s + c(-0.05, 0.05) * diff(yr_s)
  plot(sub$constant, sub$tss_slope, type = "o", col = "blue", pch = 16, lwd = 2,
       cex = 1.1, ylim = yr_s, xlab = "Constant", ylab = "Slope",
       main = paste0(ttl, "\nSlope"), cex.lab = 1.2, cex.main = 0.95, cex.axis = 1.0)
  lines(sub$constant, sub$clr_slope, type = "o", col = "red", pch = 17, lwd = 2, cex = 1.1)
  abline(h = 1, col = "grey30", lty = 2, lwd = 1.5)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)

  par(new = TRUE)
  plot(st$constant, st$median_log_tss, type = "l", col = "blue",
       lty = 3, lwd = 1.5, axes = FALSE, xlab = "", ylab = "", ylim = med_yr)
  lines(st$constant, st$median_log_clr, col = "red", lty = 3, lwd = 1.5)
  axis(4, col.axis = "grey30", cex.axis = 0.8)
  mtext(expression("Med " * -log[10](p)), side = 4, line = 2.5, cex = 0.7, col = "grey30")

  legend("topleft", legend = c("TSS side", "CLR side", "Med TSS", "Med CLR"),
         col = c("blue", "red", "blue", "red"), pch = c(16, 17, NA, NA),
         lty = c(1, 1, 3, 3), lwd = c(2, 2, 1.5, 1.5), cex = 0.65, bg = "white")

  # Right: intercept + median
  par(mar = c(5, 5, 4, 5))
  yr_i <- range(c(sub$tss_intercept, sub$clr_intercept, 0), na.rm = TRUE)
  yr_i <- yr_i + c(-0.05, 0.05) * diff(yr_i)
  plot(sub$constant, sub$tss_intercept, type = "o", col = "blue", pch = 16, lwd = 2,
       cex = 1.1, ylim = yr_i, xlab = "Constant", ylab = "Intercept",
       main = paste0(ttl, "\nIntercept"), cex.lab = 1.2, cex.main = 0.95, cex.axis = 1.0)
  lines(sub$constant, sub$clr_intercept, type = "o", col = "red", pch = 17, lwd = 2, cex = 1.1)
  abline(h = 0, col = "grey30", lty = 2, lwd = 1.5)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)

  par(new = TRUE)
  plot(st$constant, st$median_log_tss, type = "l", col = "blue",
       lty = 3, lwd = 1.5, axes = FALSE, xlab = "", ylab = "", ylim = med_yr)
  lines(st$constant, st$median_log_clr, col = "red", lty = 3, lwd = 1.5)
  axis(4, col.axis = "grey30", cex.axis = 0.8)
  mtext(expression("Med " * -log[10](p)), side = 4, line = 2.5, cex = 0.7, col = "grey30")

  par(mfrow = c(1, 1))

  # ---- Page 5: Median effect size (|estimate|) across constants ---------------
  par(mar = c(5, 5, 4, 2))
  es_yr <- range(c(st$median_es_clr, st$median_es_tss), na.rm = TRUE)
  es_yr <- es_yr + c(-0.05, 0.05) * diff(es_yr)

  plot(st$constant, st$median_es_tss, type = "o", col = "blue", pch = 16, lwd = 2,
       cex = 1.2, ylim = es_yr, xlab = "Constant",
       ylab = "Median |effect size|",
       main = paste0(ttl, "\nMedian absolute effect size across all features"),
       cex.lab = 1.3, cex.main = 1.1, cex.axis = 1.1)
  lines(st$constant, st$median_es_clr, type = "o", col = "red", pch = 17, lwd = 2, cex = 1.2)
  abline(v = 0, col = "grey50", lty = 3)
  grid(col = "grey85", lty = 1)
  legend("topleft",
         legend = c("TSS median |estimate|", "CLR median |estimate|"),
         col = c("blue", "red"), pch = c(16, 17),
         lty = 1, lwd = 2, cex = 0.9, bg = "white")

  # ---- Pages 6–26: One page per constant, regression lines -------------------
  x_seq <- seq(0, 15, length.out = 100)

  for (i in seq_len(nrow(sub))) {
    k_i     <- sub$constant[i]
    ts_s    <- sub$tss_slope[i];     ts_int <- sub$tss_intercept[i];  ts_r2 <- sub$tss_r2[i]
    cs_s    <- sub$clr_slope[i];     cs_int <- sub$clr_intercept[i];  cs_r2 <- sub$clr_r2[i]
    n_tss   <- sub$n_tss_side[i];    n_clr  <- sub$n_clr_side[i]

    y_vals <- numeric(0)
    if (!is.na(ts_s)) y_vals <- c(y_vals, ts_int + ts_s * x_seq)
    if (!is.na(cs_s)) y_vals <- c(y_vals, cs_int + cs_s * x_seq)
    if (length(y_vals) == 0) next
    yr <- range(y_vals[y_vals >= -5 & y_vals <= 25], na.rm = TRUE)
    yr <- yr + c(-0.5, 0.5)

    par(mar = c(5, 5, 5, 2))
    plot(NULL, xlim = c(0, 15), ylim = yr,
         xlab = expression(-log[10](TSS~padj)),
         ylab = expression(-log[10](CLR~padj)),
         main = paste0(ttl, "  |  constant = ", k_i),
         cex.lab = 1.3, cex.main = 1.1, cex.axis = 1.1)
    grid(col = "grey85", lty = 1)

    if (!is.na(ts_s)) lines(x_seq, ts_int + ts_s * x_seq, col = "blue", lwd = 2.5)
    if (!is.na(cs_s)) lines(x_seq, cs_int + cs_s * x_seq, col = "red", lwd = 2.5)

    tss_txt <- if (!is.na(ts_s)) {
      sprintf("TSS side (n=%d): slope=%.3f, int=%.3f, R²=%.3f", n_tss, ts_s, ts_int, ts_r2)
    } else { "TSS side: insufficient data" }

    clr_txt <- if (!is.na(cs_s)) {
      sprintf("CLR side (n=%d): slope=%.3f, int=%.3f, R²=%.3f", n_clr, cs_s, cs_int, cs_r2)
    } else { "CLR side: insufficient data" }

    legend("topleft",
           legend = c(tss_txt, clr_txt),
           col = c("blue", "red"), lty = 1, lwd = 2.5, cex = 0.85, bg = "white")
  }

  dev.off()
}

# ---- Helper: scatter plot PDF (one page per constant) -----------------------
plot_scatter_pdf <- function(pd_list, comp, label, plotdir, gamma_val) {
  comp_label <- gsub(" ", "_", comp)
  pdf(paste0(plotdir, "scatter_", label, "_", comp_label, "_gamma", gamma_val, ".pdf"),
      width = 9, height = 8)

  sig_thresh <- -log10(0.05)

  for (item in pd_list) {
    if (item$comp != comp) next
    k  <- item$k
    pd <- item$pd

    # Build title depending on which side has the constant
    if (grepl("CLR", label)) {
      ttl <- paste0(comp, " - CLR(const=", k, ") vs TSS (gamma = ", gamma_val, ")")
    } else {
      ttl <- paste0(comp, " - CLR vs TSS(const=", k, ") (gamma = ", gamma_val, ")")
    }

    n_both <- sum(pd$both);  n_tss <- sum(pd$tss_only)
    n_clr  <- sum(pd$clr_only);  n_nei <- sum(pd$neither)

    par(mar = c(5, 5, 4, 2))
    plot(pd$log_tss[pd$neither], pd$log_clr[pd$neither],
         col = "grey70", pch = 16, cex = 0.6,
         xlim = range(pd$log_tss, na.rm = TRUE),
         ylim = range(pd$log_clr, na.rm = TRUE),
         xlab = expression(-log[10](FDR)~TSS),
         ylab = expression(-log[10](FDR)~CLR),
         main = ttl, cex.lab = 1.3, cex.main = 1.0)
    points(pd$log_tss[pd$tss_only], pd$log_clr[pd$tss_only],
           col = "blue", pch = 16, cex = 0.9)
    points(pd$log_tss[pd$clr_only], pd$log_clr[pd$clr_only],
           col = "orange", pch = 16, cex = 0.9)
    points(pd$log_tss[pd$both], pd$log_clr[pd$both],
           col = "red", pch = 16, cex = 0.9)

    abline(0, 1, lty = 2, lwd = 1.5)
    abline(h = sig_thresh, col = "grey60", lty = 3, lwd = 1)
    abline(v = sig_thresh, col = "grey60", lty = 3, lwd = 1)

    legend("topleft",
           legend = c(paste0("Both (", n_both, ")"),
                      paste0("TSS only (", n_tss, ")"),
                      paste0("CLR only (", n_clr, ")"),
                      paste0("Neither (", n_nei, ")")),
           col = c("red", "blue", "orange", "grey70"),
           pch = 16, cex = 0.9, bg = "white")
  }
  dev.off()
}

# ---- Helper: effect cross plot PDF (one page per constant) ------------------
plot_effect_pdf <- function(pd_list, comp, label, plotdir, gamma_val) {
  comp_label <- gsub(" ", "_", comp)
  pdf(paste0(plotdir, "effect_", label, "_", comp_label, "_gamma", gamma_val, ".pdf"),
      width = 14, height = 7)

  for (item in pd_list) {
    if (item$comp != comp) next
    k  <- item$k
    pd <- item$pd

    # Build labels depending on which side has the constant
    if (grepl("CLR", label)) {
      clr_lab <- paste0("CLR(const=", k, ")")
      tss_lab <- "TSS"
    } else {
      clr_lab <- "CLR"
      tss_lab <- paste0("TSS(const=", k, ")")
    }

    n_both <- sum(pd$both)
    n_tss  <- sum(pd$tss_only)
    n_clr  <- sum(pd$clr_only)
    n_nei  <- sum(pd$neither)

    par(mfrow = c(1, 2))

    # ---- Left panel: CLR effect space, all 4 categories ----
    par(mar = c(5, 5, 4, 2))
    plot(pd$se_clr[pd$neither], pd$est_clr[pd$neither],
         col = "grey70", pch = 1, cex = 0.5,
         xlim = range(pd$se_clr, na.rm = TRUE),
         ylim = range(pd$est_clr, na.rm = TRUE),
         xlab = "Std Error", ylab = "Estimate",
         main = paste0(comp, " - ", clr_lab, " effect space"),
         cex.lab = 1.2, cex.main = 0.9)
    points(pd$se_clr[pd$tss_only], pd$est_clr[pd$tss_only],
           col = "blue", pch = 16, cex = 0.9)
    points(pd$se_clr[pd$clr_only], pd$est_clr[pd$clr_only],
           col = "orange", pch = 16, cex = 0.9)
    points(pd$se_clr[pd$both], pd$est_clr[pd$both],
           col = "red", pch = 16, cex = 0.9)
    abline(h = 0, lwd = 1.5)
    abline(h = c(-1, 1), lty = 3, col = "grey60")

    legend("topleft",
           legend = c(paste0("Both (", n_both, ")"),
                      paste0("TSS only (", n_tss, ")"),
                      paste0("CLR only (", n_clr, ")"),
                      paste0("Non-sig (", n_nei, ")")),
           col = c("red", "blue", "orange", "grey70"),
           pch = c(16, 16, 16, 1), cex = 0.85, bg = "white")

    # ---- Right panel: TSS effect space, all 4 categories ----
    par(mar = c(5, 5, 4, 2))
    plot(pd$se_tss[pd$neither], pd$est_tss[pd$neither],
         col = "grey70", pch = 1, cex = 0.5,
         xlim = range(pd$se_tss, na.rm = TRUE),
         ylim = range(pd$est_tss, na.rm = TRUE),
         xlab = "Std Error", ylab = "Estimate",
         main = paste0(comp, " - ", tss_lab, " effect space"),
         cex.lab = 1.2, cex.main = 0.9)
    points(pd$se_tss[pd$tss_only], pd$est_tss[pd$tss_only],
           col = "blue", pch = 16, cex = 0.9)
    points(pd$se_tss[pd$clr_only], pd$est_tss[pd$clr_only],
           col = "orange", pch = 16, cex = 0.9)
    points(pd$se_tss[pd$both], pd$est_tss[pd$both],
           col = "red", pch = 16, cex = 0.9)
    abline(h = 0, lwd = 1.5)
    abline(h = c(-1, 1), lty = 3, col = "grey60")

    legend("topleft",
           legend = c(paste0("Both (", n_both, ")"),
                      paste0("TSS only (", n_tss, ")"),
                      paste0("CLR only (", n_clr, ")"),
                      paste0("Non-sig (", n_nei, ")")),
           col = c("red", "blue", "orange", "grey70"),
           pch = c(16, 16, 16, 1), cex = 0.85, bg = "white")

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
  clr_datadir <- paste0(maindir, "CLR_constant/gamma_", gamma_val, "/data/")
  clr_plotdir <- paste0(maindir, "CLR_constant/gamma_", gamma_val, "/plots/")
  tss_datadir <- paste0(maindir, "TSS_constant/gamma_", gamma_val, "/data/")
  tss_plotdir <- paste0(maindir, "TSS_constant/gamma_", gamma_val, "/plots/")

  dir.create(clr_datadir, showWarnings = FALSE, recursive = TRUE)
  dir.create(clr_plotdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tss_datadir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tss_plotdir, showWarnings = FALSE, recursive = TRUE)

  # ###########################################################################
  # PART 1: CLR(const) vs TSS(standard)
  # ###########################################################################

  cat("\n========== PART 1: CLR constant | gamma =", gamma_val, " ==========\n")
  summary_stats1  <- data.frame()
  slope_table1    <- data.frame()
  red_dot_detail1 <- data.frame()
  plot_data_list1 <- list()

  for (k in constants) {
    cat("=== [CLR const] constant =", k, " ===\n")

    clr_const <- function(X, logComp, gamma = 0.5) {
      P  <- nrow(X); ns <- dim(logComp)[3]
      LambdaScale <- matrix(rnorm(P * ns, k, gamma), P, ns)
      logScale <- t(X) %*% LambdaScale
      return(logScale)
    }

    for (pr in pairs_list) {
      name1 <- pr[1]; name2 <- pr[2]
      comp  <- paste(name1, "vs", name2)
      cat("  ", comp, "\n")

      Y     <- cbind(datasets[[name1]], datasets[[name2]])
      conds <- c(rep(name1, ncol(datasets[[name1]])),
                 rep(name2, ncol(datasets[[name2]])))
      data  <- data.frame(condition = conds)

      res_clr <- ALDEx3::aldex(Y, ~condition, data,
                               nsample = nsample, scale = clr_const, gamma = gamma_val)
      res_tss <- ALDEx3::aldex(Y, ~condition, data,
                               nsample = nsample, scale = tss.sm, gamma = gamma_val)

      out <- run_analysis(summary(res_clr), summary(res_tss), k, comp)
      summary_stats1  <- rbind(summary_stats1,  out$stats)
      slope_table1    <- rbind(slope_table1,    out$slope)
      red_dot_detail1 <- rbind(red_dot_detail1, out$detail)
      plot_data_list1[[length(plot_data_list1) + 1]] <- list(k = k, comp = comp, pd = out$plot_data)
    }
  }

  # Save Part 1
  write.csv(summary_stats1,
            paste0(clr_datadir, "CLR_const_summary_stats_gamma", gamma_val, ".csv"),
            row.names = FALSE)
  write.csv(red_dot_detail1,
            paste0(clr_datadir, "CLR_const_red_dot_details_gamma", gamma_val, ".csv"),
            row.names = FALSE)
  save_split_rda(slope_table1, "CLR_const", clr_datadir, gamma_val)
  cat("Part 1 data saved to:", clr_datadir, "\n")

  # ###########################################################################
  # PART 2: CLR(standard) vs TSS(const)
  # ###########################################################################

  cat("\n========== PART 2: TSS constant | gamma =", gamma_val, " ==========\n")
  summary_stats2  <- data.frame()
  slope_table2    <- data.frame()
  red_dot_detail2 <- data.frame()
  plot_data_list2 <- list()

  for (k in constants) {
    cat("=== [TSS const] constant =", k, " ===\n")

    tss_const <- function(X, logComp, gamma = 0.5) {
      P  <- nrow(X); ns <- dim(logComp)[3]
      LambdaScale <- matrix(rnorm(P * ns, k, gamma), P, ns)
      logScale <- t(X) %*% LambdaScale
      return(logScale)
    }

    for (pr in pairs_list) {
      name1 <- pr[1]; name2 <- pr[2]
      comp  <- paste(name1, "vs", name2)
      cat("  ", comp, "\n")

      Y     <- cbind(datasets[[name1]], datasets[[name2]])
      conds <- c(rep(name1, ncol(datasets[[name1]])),
                 rep(name2, ncol(datasets[[name2]])))
      data  <- data.frame(condition = conds)

      res_clr <- ALDEx3::aldex(Y, ~condition, data,
                               nsample = nsample, scale = clr.sm, gamma = gamma_val)
      res_tss <- ALDEx3::aldex(Y, ~condition, data,
                               nsample = nsample, scale = tss_const, gamma = gamma_val)

      out <- run_analysis(summary(res_clr), summary(res_tss), k, comp)
      summary_stats2  <- rbind(summary_stats2,  out$stats)
      slope_table2    <- rbind(slope_table2,    out$slope)
      red_dot_detail2 <- rbind(red_dot_detail2, out$detail)
      plot_data_list2[[length(plot_data_list2) + 1]] <- list(k = k, comp = comp, pd = out$plot_data)
    }
  }

  # Save Part 2
  write.csv(summary_stats2,
            paste0(tss_datadir, "TSS_const_summary_stats_gamma", gamma_val, ".csv"),
            row.names = FALSE)
  write.csv(red_dot_detail2,
            paste0(tss_datadir, "TSS_const_red_dot_details_gamma", gamma_val, ".csv"),
            row.names = FALSE)
  save_split_rda(slope_table2, "TSS_const", tss_datadir, gamma_val)
  cat("Part 2 data saved to:", tss_datadir, "\n")

  # ###########################################################################
  # PLOTS
  # ###########################################################################

  for (comp in unique(slope_table1$comparison)) {
    plot_pairwise(slope_table1, summary_stats1, comp, "CLR_const", clr_plotdir, gamma_val)
    plot_scatter_pdf(plot_data_list1, comp, "CLR_const", clr_plotdir, gamma_val)
    plot_effect_pdf(plot_data_list1, comp, "CLR_const", clr_plotdir, gamma_val)
  }
  for (comp in unique(slope_table2$comparison)) {
    plot_pairwise(slope_table2, summary_stats2, comp, "TSS_const", tss_plotdir, gamma_val)
    plot_scatter_pdf(plot_data_list2, comp, "TSS_const", tss_plotdir, gamma_val)
    plot_effect_pdf(plot_data_list2, comp, "TSS_const", tss_plotdir, gamma_val)
  }

  cat("Plots saved for gamma =", gamma_val, "\n")
}

cat("\n========== All done! ==========\n")
cat("Main folder:", maindir, "\n\n")
cat("Structure: {CLR,TSS}_constant/gamma_{0.1,0.3,0.5}/{data,plots}/\n")
