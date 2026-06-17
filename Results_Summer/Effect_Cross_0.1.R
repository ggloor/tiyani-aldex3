# effect_cross_aldex_edgeR.R
# Effect plots with cross-method significance overlay.
#
# ALDEx3 side: Std Error (x) vs Estimate (y) on ALL 1117 OTUs,
#              with edgeR significance overlaid for the filterByExpr subset.
# edgeR side:  logCPM (x) vs logFC (y) on ALL 1117 OTUs,
#              unfiltered edgeR run for coordinates, filtered FDR for significance.
#
# Uses gamma = 0.3. Boundary lines shown as light grey dotted lines.

library(ALDEx3)
library(edgeR)
set.seed(12345)

# --- Paths and parameters ---
data_dir <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/data"
out_dir  <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer/aldex_vs_edger/effect_cross_edgeR0.1"
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

nsample   <- 128
gamma_val <- 0.3
fdr_cut   <- 0.05

# --- Load cohorts ---
load(file.path(data_dir, "kin.Rda"))
load(file.path(data_dir, "pup.Rda"))
load(file.path(data_dir, "mid.Rda"))
load(file.path(data_dir, "you.Rda"))
load(file.path(data_dir, "mage.Rda"))
load(file.path(data_dir, "eld.Rda"))
load(file.path(data_dir, "cent.Rda"))

datasets <- list(kin = kin, pup = pup, mid = mid, you = you,
                 mage = mage, eld = eld, cent = cent)
pairs <- combn(names(datasets), 2)


# #############################################################################
# PART 1: DATA ANALYSIS
# #############################################################################

cat("===== PART 1: ANALYSIS =====\n\n")

all_results <- list()

for (i in 1:ncol(pairs)) {
  
  name1 <- pairs[1, i]
  name2 <- pairs[2, i]
  comp  <- paste(name1, "vs", name2)
  cat(sprintf("[%2d/21] %s\n", i, comp))
  
  Y     <- cbind(datasets[[name1]], datasets[[name2]])
  conds <- c(rep(name1, ncol(datasets[[name1]])),
             rep(name2, ncol(datasets[[name2]])))
  dd    <- data.frame(condition = conds)
  group <- factor(conds)
  
  # --- ALDEx3 CLR (full OTU set) ---
  cat("  ALDEx3 CLR\n")
  res_clr <- ALDEx3::aldex(Y, ~condition, dd,
                           nsample = nsample, scale = clr.sm, gamma = gamma_val)
  sum_clr <- summary(res_clr)
  
  # --- ALDEx3 TSS (full OTU set) ---
  cat("  ALDEx3 TSS\n")
  res_tss <- ALDEx3::aldex(Y, ~condition, dd,
                           nsample = nsample, scale = tss.sm, gamma = gamma_val)
  sum_tss <- summary(res_tss)
  
  # --- edgeR UNFILTERED (for logCPM/logFC coordinates on all OTUs) ---
  cat("  edgeR TMM unfiltered (for coordinates)\n")
  dge_all    <- DGEList(counts = Y, group = group)
  dge_all    <- normLibSizes(dge_all)
  design     <- model.matrix(~group)
  fit_all    <- glmQLFit(dge_all, design)
  qlf_all    <- glmQLFTest(fit_all, coef = 2)
  tt_all     <- topTags(qlf_all, n = Inf, sort.by = "none")$table
  
  # compute per-OTU standard error from the QL fit
  edger_se_all <- abs(tt_all$logFC) / sqrt(pmax(tt_all$F, 1e-10))
  
  # --- edgeR FILTERED (for significance calls) ---
  cat("  edgeR TMM filtered (for FDR)\n")
  dge_f  <- DGEList(counts = Y, group = group)
  keep   <- filterByExpr(dge_f, group = group)
  cat(sprintf("  filterByExpr kept %d / %d OTUs\n", sum(keep), nrow(Y)))
  
  dge_f  <- dge_f[keep, , keep.lib.sizes = FALSE]
  dge_f  <- normLibSizes(dge_f)
  fit_f  <- glmQLFit(dge_f, design)
  qlf_f  <- glmQLFTest(fit_f, coef = 2)
  tt_f   <- topTags(qlf_f, n = Inf, sort.by = "none")$table
  
  # --- Build filtered FDR expanded to all 1117 OTUs ---
  edger_fdr_full <- rep(1, nrow(Y))
  names(edger_fdr_full) <- rownames(Y)
  edger_fdr_full[keep] <- tt_f$FDR
  
  # --- Save ---
  result <- list(
    # Full ALDEx3 results (all 1117 OTUs)
    clr_full = data.frame(OTU       = rownames(sum_clr),
                          estimate  = sum_clr$estimate,
                          std_error = sum_clr$std.error,
                          pval_adj  = sum_clr$p.val.adj),
    tss_full = data.frame(OTU       = rownames(sum_tss),
                          estimate  = sum_tss$estimate,
                          std_error = sum_tss$std.error,
                          pval_adj  = sum_tss$p.val.adj),
    # edgeR on ALL OTUs (for effect plot coordinates)
    edger_all = data.frame(OTU       = rownames(tt_all),
                           logFC     = tt_all$logFC,
                           logCPM    = tt_all$logCPM,
                           std_error = edger_se_all),
    # Filtered FDR expanded to all 1117 (filtered-out = 1)
    edger_fdr_full = edger_fdr_full,
    comparison = comp,
    n_kept     = sum(keep),
    keep       = keep
  )
  all_results[[comp]] <- result
  
  fname <- paste0(name1, "_vs_", name2, ".Rda")
  save(result, file = file.path(out_dir, fname))
}

cat("\nAll results saved to:", out_dir, "\n\n")


# #############################################################################
# PART 2: PLOTTING
# #############################################################################

cat("===== PART 2: PLOTS =====\n\n")

# --- ALDEx3 effect plot (all 1117 OTUs, overlay = other method's sig only) ---
plot_effect_aldex <- function(std_error, estimate, own_padj, other_fdr,
                              own_label, other_label, comp_name) {
  
  sig_own   <- own_padj  < fdr_cut
  sig_other <- other_fdr < fdr_cut
  
  n_both       <- sum(sig_own & sig_other)
  n_other_only <- sum(!sig_own & sig_other)
  n_nonsig     <- sum(!sig_other)   # everything edgeR doesn't call
  
  # only color the OTHER method's significant calls
  cols <- rep("grey80", length(estimate))
  cols[sig_own & sig_other]  <- "red"      # both agree
  cols[!sig_own & sig_other] <- "orange"   # other method only
  
  # plot non-sig first, sig on top
  ord <- order(cols != "grey80")
  
  plot(std_error[ord], estimate[ord],
       xlab = "Std Error",
       ylab = "Estimate",
       main = sprintf("%s \u2013 %s (%s sig overlay)", comp_name, own_label, other_label),
       pch = 19, cex = 0.5, col = cols[ord])
  
  abline(h = 0, col = "black", lwd = 1)
  
  # boundary lines
  abline(h =  1, col = "grey70", lty = 3, lwd = 0.8)
  abline(h = -1, col = "grey70", lty = 3, lwd = 0.8)
  
  legend("topleft",
         legend = c(sprintf("Both (%d)", n_both),
                    sprintf("%s only (%d)", other_label, n_other_only),
                    sprintf("Non-sig (%d)", n_nonsig)),
         col = c("red", "orange", "grey80"),
         pch = 19, cex = 0.75, bty = "n")
}


# --- edgeR effect plot (all 1117 OTUs, overlay = other method's sig only) ---
plot_effect_edger <- function(std_error, logFC, edger_fdr_full, other_padj,
                              own_label, other_label, comp_name) {
  
  sig_own   <- edger_fdr_full < fdr_cut
  sig_other <- other_padj     < fdr_cut
  
  n_both       <- sum(sig_own & sig_other)
  n_other_only <- sum(!sig_own & sig_other)
  n_nonsig     <- sum(!sig_other)   # everything ALDEx3 doesn't call
  
  # only color the OTHER method's significant calls
  cols <- rep("grey80", length(logFC))
  cols[sig_own & sig_other]  <- "red"      # both agree
  cols[!sig_own & sig_other] <- "orange"   # other method only
  
  ord <- order(cols != "grey80")
  
  plot(std_error[ord], logFC[ord],
       xlab = "Std Error",
       ylab = "logFC",
       main = sprintf("%s \u2013 %s (%s sig overlay)", comp_name, own_label, other_label),
       pch = 19, cex = 0.5, col = cols[ord])
  
  abline(h = 0, col = "black", lwd = 1)
  
  # boundary lines
  abline(h =  1, col = "grey70", lty = 3, lwd = 0.8)
  abline(h = -1, col = "grey70", lty = 3, lwd = 0.8)
  
  legend("topleft",
         legend = c(sprintf("Both (%d)", n_both),
                    sprintf("%s only (%d)", other_label, n_other_only),
                    sprintf("Non-sig (%d)", n_nonsig)),
         col = c("red", "orange", "grey80"),
         pch = 19, cex = 0.75, bty = "n")
}


# =========================================================================
# PLOT SET 1: CLR <-> edgeR
# =========================================================================
cat("--- Plot Set 1: CLR <-> edgeR ---\n")

pdf(file.path(out_dir, "plots", "effect_cross_CLR_edgeR_gamma03.pdf"),
    width = 10, height = 5)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- all_results[[comp]]
  
  par(mfrow = c(1, 2))
  
  # left: CLR effect plot (all 1117), edgeR sig overlay
  plot_effect_aldex(r$clr_full$std_error, r$clr_full$estimate,
                    r$clr_full$pval_adj, r$edger_fdr_full,
                    "CLR", "edgeR TMM", comp)
  
  # right: edgeR effect plot (all 1117), CLR sig overlay
  plot_effect_edger(r$edger_all$std_error, r$edger_all$logFC,
                    r$edger_fdr_full, r$clr_full$pval_adj,
                    "edgeR TMM", "CLR", comp)
}
dev.off()
cat("  Saved: effect_cross_CLR_edgeR_gamma03.pdf\n")


# =========================================================================
# PLOT SET 2: TSS <-> edgeR
# =========================================================================
cat("--- Plot Set 2: TSS <-> edgeR ---\n")

pdf(file.path(out_dir, "plots", "effect_cross_TSS_edgeR_gamma03.pdf"),
    width = 10, height = 5)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- all_results[[comp]]
  
  par(mfrow = c(1, 2))
  
  # left: TSS effect plot (all 1117), edgeR sig overlay
  plot_effect_aldex(r$tss_full$std_error, r$tss_full$estimate,
                    r$tss_full$pval_adj, r$edger_fdr_full,
                    "TSS", "edgeR TMM", comp)
  
  # right: edgeR effect plot (all 1117), TSS sig overlay
  plot_effect_edger(r$edger_all$std_error, r$edger_all$logFC,
                    r$edger_fdr_full, r$tss_full$pval_adj,
                    "edgeR TMM", "TSS", comp)
}
dev.off()
cat("  Saved: effect_cross_TSS_edgeR_gamma03.pdf\n")


# =========================================================================
# PLOT SET 3: All four panels
# =========================================================================
cat("--- Plot Set 3: all four panels ---\n")

pdf(file.path(out_dir, "plots", "effect_cross_all_panels_gamma03.pdf"),
    width = 10, height = 10)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- all_results[[comp]]
  
  par(mfrow = c(2, 2))
  
  plot_effect_aldex(r$clr_full$std_error, r$clr_full$estimate,
                    r$clr_full$pval_adj, r$edger_fdr_full,
                    "CLR", "edgeR TMM", comp)
  
  plot_effect_edger(r$edger_all$std_error, r$edger_all$logFC,
                    r$edger_fdr_full, r$clr_full$pval_adj,
                    "edgeR TMM", "CLR", comp)
  
  plot_effect_aldex(r$tss_full$std_error, r$tss_full$estimate,
                    r$tss_full$pval_adj, r$edger_fdr_full,
                    "TSS", "edgeR TMM", comp)
  
  plot_effect_edger(r$edger_all$std_error, r$edger_all$logFC,
                    r$edger_fdr_full, r$tss_full$pval_adj,
                    "edgeR TMM", "TSS", comp)
}
dev.off()
cat("  Saved: effect_cross_all_panels_gamma03.pdf\n")


# #############################################################################
# PART 3: SUMMARY TABLE
# #############################################################################

cat("\n===== PART 3: SUMMARY =====\n\n")

cat(sprintf("=== Cross-method agreement (FDR < 0.05, gamma = %g) ===\n", gamma_val))
cat("  ALDEx3 on all 1117 OTUs; edgeR significance from filtered subset\n\n")
cat(sprintf("%-18s %6s %6s %6s %6s %8s %8s %9s %9s\n",
            "Comparison", "Kept", "CLR", "TSS", "edgeR",
            "CLR+edgR", "TSS+edgR", "edgR only", "edgR only"))
cat(sprintf("%-18s %6s %6s %6s %6s %8s %8s %9s %9s\n",
            "", "", "(all)", "(all)", "(filt)",
            "", "", "(vs CLR)", "(vs TSS)"))
cat(paste(rep("-", 88), collapse = ""), "\n")

for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- all_results[[comp]]
  
  sig_clr   <- r$clr_full$pval_adj < fdr_cut
  sig_tss   <- r$tss_full$pval_adj < fdr_cut
  sig_edger <- r$edger_fdr_full    < fdr_cut
  
  cat(sprintf("%-18s %6d %6d %6d %6d %8d %8d %9d %9d\n", comp, r$n_kept,
              sum(sig_clr), sum(sig_tss), sum(sig_edger),
              sum(sig_clr & sig_edger),
              sum(sig_tss & sig_edger),
              sum(sig_edger & !sig_clr),
              sum(sig_edger & !sig_tss)))
}

cat("\nDone.\n")