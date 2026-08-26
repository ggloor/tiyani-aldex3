# =============================================================================
# TEST: coefficient.sm effect space — single run, no replicates
#
# Runs coefficient.sm(c.mu=c(0,k)) vs CLR(std) vs TSS for one constant (k=0.3)
# and one gamma (0.3). Produces:
#   Plot 1: 3-panel effect cross (unified 5 categories across all panels)
#   Plot 2: 2-panel -log10(FDR) scatter (Coeff vs CLR, Coeff vs TSS)
#
# Uses coefficient.sm with c.sd (std dev per coefficient) per ALDEx3 v1.3.0 API:
#   c.mu = c(intercept_mean, slope_mean)  → c(0, k)
#   c.sd = c(intercept_sd, slope_sd)      → c(0, gamma)
#
# Robustness = all_three / (all_three + coeff_only + clr_only + tss_only)
#   1 = all normalizations agree, 0 = max disagreement
#
# Zero p-value handling: min(nonzero)/10 before -log10 (Dr. Gloor convention)
# =============================================================================

library(ALDEx3)

# ---- Config -----------------------------------------------------------------
datadir   <- "~/Desktop/3383/0_git/tiyani-aldex3/data/"
nsample   <- 32
gamma_val <- 0.3
k         <- 0.3

# ---- Output -----------------------------------------------------------------
outdir <- paste0("/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/coefficient_sm_TEST/")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data --------------------------------------------------------------
load(paste0(datadir, "cent.Rda"))
load(paste0(datadir, "mage.Rda"))

Y     <- cbind(cent, mage)
conds <- c(rep("cent", ncol(cent)), rep("mage", ncol(mage)))
dat   <- data.frame(condition = conds)

# ---- Run all three methods --------------------------------------------------
cat("Running coefficient.sm (c.mu=c(0,", k, "), c.sd=c(0,", gamma_val, "))...\n")
res_coeff <- ALDEx3::aldex(Y, ~condition, dat,
                           nsample = nsample, scale = coefficient.sm,
                           c.mu = c(0, k),
                           c.sd = c(0, gamma_val))

cat("Running CLR (standard)...\n")
res_clr <- ALDEx3::aldex(Y, ~condition, dat,
                         nsample = nsample, scale = clr.sm, gamma = gamma_val)

cat("Running TSS (standard)...\n")
res_tss <- ALDEx3::aldex(Y, ~condition, dat,
                         nsample = nsample, scale = tss.sm, gamma = gamma_val)

sum_coeff <- summary(res_coeff)
sum_clr   <- summary(res_clr)
sum_tss   <- summary(res_tss)

# ---- Significance -----------------------------------------------------------
sig_coeff <- sum_coeff$p.val.adj < 0.05
sig_clr   <- sum_clr$p.val.adj < 0.05
sig_tss   <- sum_tss$p.val.adj < 0.05

# ---- 5 unified categories (same logic as effect sweep & 3.0 scripts) -------
# Priority: all3 > coeff > clr > tss > none
n    <- nrow(sum_coeff)
cats <- rep("none", n)
cats[sig_tss & !sig_clr & !sig_coeff]   <- "tss"
cats[sig_clr & !sig_tss & !sig_coeff]   <- "clr"
cats[sig_coeff & !sig_clr & !sig_tss]   <- "coeff"
cats[sig_clr & sig_tss & !sig_coeff]    <- "clr"    # CLR+TSS agree, coeff doesn't
cats[sig_coeff & sig_clr & !sig_tss]    <- "coeff"  # coeff+CLR, no TSS
cats[sig_coeff & sig_tss & !sig_clr]    <- "coeff"  # coeff+TSS, no CLR
cats[sig_coeff & sig_clr & sig_tss]     <- "all3"

n_all3  <- sum(cats == "all3")
n_coeff <- sum(cats == "coeff")
n_clr   <- sum(cats == "clr")
n_tss   <- sum(cats == "tss")
n_none  <- sum(cats == "none")

# Robustness: 1 = all methods agree, 0 = max disagreement
n_any_sig  <- n_all3 + n_coeff + n_clr + n_tss
robustness <- if (n_any_sig > 0) round(n_all3 / n_any_sig, 3) else NA

cat("\n--- Counts ---\n")
cat("Coeff sig:", sum(sig_coeff), "\n")
cat("CLR sig:  ", sum(sig_clr), "\n")
cat("TSS sig:  ", sum(sig_tss), "\n")
cat("All three:", n_all3, "\n")
cat("Coeff only:", n_coeff, "\n")
cat("CLR only:  ", n_clr, "\n")
cat("TSS only:  ", n_tss, "\n")
cat("Non-sig:   ", n_none, "\n")
cat("Robustness:", robustness, "\n")

# ---- Helper: plot one effect panel with all 5 categories --------------------
plot_effect_panel <- function(sum_df, cats, title_text, method_sig,
                              cnt_all3, cnt_coeff, cnt_clr, cnt_tss, cnt_none) {

  # Plot non-sig first (background)
  plot(sum_df$std.error[cats == "none"], sum_df$estimate[cats == "none"],
       col = "grey75", pch = 1, cex = 0.4,
       xlim = range(sum_df$std.error), ylim = range(sum_df$estimate),
       xlab = "Std Error", ylab = "Effect",
       main = title_text, cex.lab = 1.2, cex.main = 0.9)

  # Overlay each category
  if (any(cats == "tss"))
    points(sum_df$std.error[cats == "tss"], sum_df$estimate[cats == "tss"],
           col = "forestgreen", pch = 16, cex = 0.8)
  if (any(cats == "clr"))
    points(sum_df$std.error[cats == "clr"], sum_df$estimate[cats == "clr"],
           col = "dodgerblue", pch = 16, cex = 0.8)
  if (any(cats == "coeff"))
    points(sum_df$std.error[cats == "coeff"], sum_df$estimate[cats == "coeff"],
           col = "orange", pch = 16, cex = 0.8)
  if (any(cats == "all3"))
    points(sum_df$std.error[cats == "all3"], sum_df$estimate[cats == "all3"],
           col = "red", pch = 16, cex = 0.8)

  abline(h = 0, lwd = 1.5)
  abline(h = c(-1, 1), lty = 3, col = "grey60")

  # Per-panel legend with panel-specific sig count
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

# ---- Generate plots (saved to PDF) ------------------------------------------
pdf_path <- paste0(outdir, "coeff_test_k", k, "_gamma", gamma_val, ".pdf")
pdf(pdf_path, width = 16, height = 7)

# ---- Plot 1: 3-panel effect cross (Coeff, CLR, TSS) — unified categories ---
par(mfrow = c(1, 3), mar = c(5, 5, 3, 1), oma = c(0, 0, 3, 0))

plot_effect_panel(sum_coeff, cats, "Coefficient.sm effect space", sum(sig_coeff),
                  n_all3, n_coeff, n_clr, n_tss, n_none)
plot_effect_panel(sum_clr, cats, "CLR(std) effect space", sum(sig_clr),
                  n_all3, n_coeff, n_clr, n_tss, n_none)
plot_effect_panel(sum_tss, cats, "TSS(std) effect space", sum(sig_tss),
                  n_all3, n_coeff, n_clr, n_tss, n_none)

# Overall title
mtext(paste0("cent vs mage  |  k = ", k, "   gamma = ", gamma_val,
             "   |   All three: ", n_all3,
             "   Robustness: ", robustness),
      side = 3, outer = TRUE, cex = 1.0, font = 2)

# ---- Plot 2: 2-panel scatter (-log10 FDR) ----------------------------------

# Zero p-value handling: min(nonzero)/10
coeff_padj <- sum_coeff$p.val.adj
clr_padj   <- sum_clr$p.val.adj
tss_padj   <- sum_tss$p.val.adj
all_pvals  <- c(coeff_padj, clr_padj, tss_padj)
min_nz     <- min(all_pvals[all_pvals > 0])
coeff_padj[coeff_padj == 0] <- min_nz / 10
clr_padj[clr_padj == 0]     <- min_nz / 10
tss_padj[tss_padj == 0]     <- min_nz / 10

log_coeff  <- -log10(coeff_padj)
log_clr    <- -log10(clr_padj)
log_tss    <- -log10(tss_padj)
sig_thresh <- -log10(0.05)

# Color vector matching cats
scatter_cols <- rep("grey75", n)
scatter_cols[cats == "tss"]   <- "forestgreen"
scatter_cols[cats == "clr"]   <- "dodgerblue"
scatter_cols[cats == "coeff"] <- "orange"
scatter_cols[cats == "all3"]  <- "red"
scatter_pch <- ifelse(cats == "none", 1, 16)
scatter_cex <- ifelse(cats == "none", 0.4, 0.7)

par(mfrow = c(1, 2), mar = c(5, 5, 3, 2), oma = c(0, 0, 3, 0))

# Left: Coeff vs CLR scatter
xy_range <- range(c(log_clr, log_coeff))
plot(log_clr, log_coeff,
     col = scatter_cols, pch = scatter_pch, cex = scatter_cex,
     xlim = xy_range, ylim = xy_range,
     xlab = expression(-log[10](FDR) ~ CLR),
     ylab = expression(-log[10](FDR) ~ Coefficient),
     main = "Coeff vs CLR",
     cex.lab = 1.2, cex.main = 0.95)
abline(0, 1, lty = 2, lwd = 1.5)
abline(h = sig_thresh, col = "grey60", lty = 3)
abline(v = sig_thresh, col = "grey60", lty = 3)
legend("topleft",
       legend = c(paste0("All three (", n_all3, ")"),
                  paste0("Coeff only (", n_coeff, ")"),
                  paste0("CLR only (", n_clr, ")"),
                  paste0("TSS only (", n_tss, ")"),
                  paste0("Non-sig (", n_none, ")")),
       col = c("red", "orange", "dodgerblue", "forestgreen", "grey75"),
       pch = c(16, 16, 16, 16, 1), cex = 0.75, bg = "white")

# Right: Coeff vs TSS scatter
xy_range2 <- range(c(log_tss, log_coeff))
plot(log_tss, log_coeff,
     col = scatter_cols, pch = scatter_pch, cex = scatter_cex,
     xlim = xy_range2, ylim = xy_range2,
     xlab = expression(-log[10](FDR) ~ TSS),
     ylab = expression(-log[10](FDR) ~ Coefficient),
     main = "Coeff vs TSS",
     cex.lab = 1.2, cex.main = 0.95)
abline(0, 1, lty = 2, lwd = 1.5)
abline(h = sig_thresh, col = "grey60", lty = 3)
abline(v = sig_thresh, col = "grey60", lty = 3)
legend("topleft",
       legend = c(paste0("All three (", n_all3, ")"),
                  paste0("Coeff only (", n_coeff, ")"),
                  paste0("CLR only (", n_clr, ")"),
                  paste0("TSS only (", n_tss, ")"),
                  paste0("Non-sig (", n_none, ")")),
       col = c("red", "orange", "dodgerblue", "forestgreen", "grey75"),
       pch = c(16, 16, 16, 16, 1), cex = 0.75, bg = "white")

# Overall title
mtext(paste0("cent vs mage  |  k = ", k, "   gamma = ", gamma_val,
             "   |   Robustness: ", robustness),
      side = 3, outer = TRUE, cex = 1.0, font = 2)

par(mfrow = c(1, 1))

dev.off()
cat("Plots saved to:", pdf_path, "\n")
cat("\n========== TEST done! ==========\n")
