# Plot_Filtered_Coords.R
# ============================================================================
# Effect cross plots: ALDEx3 (CLR/TSS) vs edgeR (TMM, TMMwsp, RLE)
# using FILTERED edgeR coordinates (logFC/logCPM from filterByExpr-kept OTUs).
#
# This script differs from Updated_Results0.1.R PART 2 in one key way:
#   - The edgeR MA panels use coords_filt (filtered run) instead of
#     coords_unfilt (unfiltered run).
#   - Only the kept OTUs appear on the edgeR panel.
#   - ALDEx3 panels are unchanged (always use all 1117 OTUs).
#   - Significance masks are subset to kept OTUs on the edgeR panel.
#
# Requires: .Rda files produced by Updated_Results0.1.R (with both
#           coords_unfilt and coords_filt stored).
# ============================================================================


# #############################################################################
# SECTION 1: SETUP
# #############################################################################

rda_dir  <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/aldex_vs_edger/effect_cross_edgeR_allnorms"
plot_dir <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/aldex_vs_edger/effect_cross_edgeR_filtered_coords/plots0.1"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cut      <- 0.05
norm_methods <- c("TMM", "TMMwsp", "RLE")
gamma_val    <- 0.3

# Cohort names and all 21 pairwise comparisons
cohorts <- c("kin", "pup", "mid", "you", "mage", "eld", "cent")
pairs   <- combn(cohorts, 2)


# #############################################################################
# SECTION 2: LOAD ALL .Rda FILES
# #############################################################################

cat("===== Loading .Rda files =====\n\n")

all_results <- list()

for (i in 1:ncol(pairs)) {
  name1 <- pairs[1, i]
  name2 <- pairs[2, i]
  comp  <- paste(name1, "vs", name2)
  fname <- file.path(rda_dir, paste0(name1, "_vs_", name2, ".Rda"))

  if (!file.exists(fname)) {
    warning(sprintf("Missing: %s — skipping\n", fname))
    next
  }

  load(fname)  # loads 'result'
  all_results[[comp]] <- result
  cat(sprintf("  Loaded: %s\n", comp))
}

cat(sprintf("\nLoaded %d comparisons.\n\n", length(all_results)))


# #############################################################################
# SECTION 3: CLR vs edgeR — FILTERED edgeR COORDINATES
# #############################################################################
#
# Left panel:  ALDEx3 CLR effect plot (all 1117 OTUs), edgeR sig overlay
#              — IDENTICAL to unfiltered version
# Right panel: edgeR MA plot using FILTERED coords (kept OTUs only),
#              CLR sig overlay (also subset to kept OTUs)
# #############################################################################

cat("===== Plotting CLR vs edgeR (filtered coords) =====\n\n")

for (nm in norm_methods) {

  cat(sprintf("--- CLR vs edgeR %s (filtered) ---\n", nm))

  pdf(file.path(plot_dir,
                sprintf("effect_cross_CLR_edgeR_%s_filtered_gamma%s.pdf",
                        nm, gsub("\\.", "", as.character(gamma_val)))),
      width = 12, height = 6)

  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- all_results[[comp]]

    # ---- Significance masks (all 1117 OTUs) ----
    sig_clr   <- r$clr$pval_adj    < fdr_cut
    sig_edger <- r$edger[[nm]]$fdr < fdr_cut

    both       <- sig_clr & sig_edger
    edger_only <- sig_edger & !sig_clr
    clr_only   <- sig_clr & !sig_edger
    neither    <- !sig_clr & !sig_edger

    # ---- Left panel: CLR effect plot (all OTUs, unchanged) ----
    fin_clr  <- is.finite(r$clr$estimate) & is.finite(r$clr$std_error)
    xlim_clr <- range(r$clr$std_error[fin_clr], na.rm = TRUE) * c(0.95, 1.05)
    ylim_clr <- range(r$clr$estimate[fin_clr],  na.rm = TRUE) * 1.05

    # ---- Right panel: edgeR MA using FILTERED coords (kept OTUs only) ----
    keep       <- r$edger[[nm]]$keep       # logical mask, length 1117
    filt_fc    <- r$edger[[nm]]$coords_filt$logFC
    filt_cpm   <- r$edger[[nm]]$coords_filt$logCPM
    fin_edg_f  <- is.finite(filt_fc) & is.finite(filt_cpm)

    # Subset significance masks to kept OTUs only
    both_k       <- both[keep]
    edger_only_k <- edger_only[keep]
    clr_only_k   <- clr_only[keep]
    neither_k    <- neither[keep]

    xlim_edg <- range(filt_cpm[fin_edg_f], na.rm = TRUE) * c(0.95, 1.05)
    ylim_edg <- range(filt_fc[fin_edg_f],  na.rm = TRUE) * 1.05

    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

    # ---- LEFT: CLR effect plot, edgeR sig overlay (all OTUs) ----
    plot(r$clr$std_error[neither & fin_clr], r$clr$estimate[neither & fin_clr],
         xlab = "Std Error", ylab = "Estimate",
         main = paste0(comp, " – CLR (edgeR ", nm, " sig overlay)"),
         pch = 1, cex = 0.4, col = rgb(0.5, 0.5, 0.5, 0.2),
         xlim = xlim_clr, ylim = ylim_clr)
    abline(h = 0, col = "black")
    abline(h = c(-1, 1), col = "grey60", lty = 3)

    if (any(both & fin_clr))
      points(r$clr$std_error[both & fin_clr], r$clr$estimate[both & fin_clr],
             pch = 19, cex = 0.7, col = "red")
    if (any(edger_only & fin_clr))
      points(r$clr$std_error[edger_only & fin_clr], r$clr$estimate[edger_only & fin_clr],
             pch = 19, cex = 0.7, col = "blue")

    legend("topleft",
           legend = c(paste0("Both (", sum(both), ")"),
                      paste0("edgeR only (", sum(edger_only), ")"),
                      paste0("Non-sig (", sum(!sig_edger), ")")),
           col = c("red", "blue", rgb(0.5, 0.5, 0.5, 0.4)),
           pch = c(19, 19, 1), cex = 0.7, bg = "white")

    # ---- RIGHT: edgeR MA plot — FILTERED coords, CLR sig overlay ----
    plot(filt_cpm[neither_k & fin_edg_f],
         filt_fc[neither_k & fin_edg_f],
         xlab = "logCPM (filtered)", ylab = "logFC (filtered)",
         main = paste0(comp, " – edgeR ", nm,
                       " FILTERED (", sum(keep), " OTUs, CLR sig overlay)"),
         pch = 1, cex = 0.4, col = rgb(0.5, 0.5, 0.5, 0.2),
         xlim = xlim_edg, ylim = ylim_edg)
    abline(h = 0, col = "black")
    abline(h = c(-1, 1), col = "grey60", lty = 3)

    if (any(both_k & fin_edg_f))
      points(filt_cpm[both_k & fin_edg_f],
             filt_fc[both_k & fin_edg_f],
             pch = 19, cex = 0.7, col = "red")
    if (any(clr_only_k & fin_edg_f))
      points(filt_cpm[clr_only_k & fin_edg_f],
             filt_fc[clr_only_k & fin_edg_f],
             pch = 19, cex = 0.7, col = "orange")

    legend("topleft",
           legend = c(paste0("Both (", sum(both_k), ")"),
                      paste0("CLR only (", sum(clr_only_k), ")"),
                      paste0("Non-sig (", sum(neither_k), ")")),
           col = c("red", "orange", rgb(0.5, 0.5, 0.5, 0.4)),
           pch = c(19, 19, 1), cex = 0.7, bg = "white")

    par(mfrow = c(1, 1))
  }
  dev.off()
  cat(sprintf("  Saved: CLR vs edgeR %s filtered effect cross plots\n", nm))
}


# #############################################################################
# SECTION 4: TSS vs edgeR — FILTERED edgeR COORDINATES
# #############################################################################
#
# Left panel:  ALDEx3 TSS effect plot (all 1117 OTUs), edgeR sig overlay
#              — IDENTICAL to unfiltered version
# Right panel: edgeR MA plot using FILTERED coords (kept OTUs only),
#              TSS sig overlay (also subset to kept OTUs)
# #############################################################################

cat("\n===== Plotting TSS vs edgeR (filtered coords) =====\n\n")

for (nm in norm_methods) {

  cat(sprintf("--- TSS vs edgeR %s (filtered) ---\n", nm))

  pdf(file.path(plot_dir,
                sprintf("effect_cross_TSS_edgeR_%s_filtered_gamma%s.pdf",
                        nm, gsub("\\.", "", as.character(gamma_val)))),
      width = 12, height = 6)

  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- all_results[[comp]]

    # ---- Significance masks (all 1117 OTUs) ----
    sig_tss   <- r$tss$pval_adj    < fdr_cut
    sig_edger <- r$edger[[nm]]$fdr < fdr_cut

    both       <- sig_tss & sig_edger
    edger_only <- sig_edger & !sig_tss
    tss_only   <- sig_tss & !sig_edger
    neither    <- !sig_tss & !sig_edger

    # ---- Left panel: TSS effect plot (all OTUs, unchanged) ----
    fin_tss  <- is.finite(r$tss$estimate) & is.finite(r$tss$std_error)
    xlim_tss <- range(r$tss$std_error[fin_tss], na.rm = TRUE) * c(0.95, 1.05)
    ylim_tss <- range(r$tss$estimate[fin_tss],  na.rm = TRUE) * 1.05

    # ---- Right panel: edgeR MA using FILTERED coords (kept OTUs only) ----
    keep       <- r$edger[[nm]]$keep
    filt_fc    <- r$edger[[nm]]$coords_filt$logFC
    filt_cpm   <- r$edger[[nm]]$coords_filt$logCPM
    fin_edg_f  <- is.finite(filt_fc) & is.finite(filt_cpm)

    # Subset significance masks to kept OTUs only
    both_k       <- both[keep]
    edger_only_k <- edger_only[keep]
    tss_only_k   <- tss_only[keep]
    neither_k    <- neither[keep]

    xlim_edg <- range(filt_cpm[fin_edg_f], na.rm = TRUE) * c(0.95, 1.05)
    ylim_edg <- range(filt_fc[fin_edg_f],  na.rm = TRUE) * 1.05

    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

    # ---- LEFT: TSS effect plot, edgeR sig overlay (all OTUs) ----
    plot(r$tss$std_error[neither & fin_tss], r$tss$estimate[neither & fin_tss],
         xlab = "Std Error", ylab = "Estimate",
         main = paste0(comp, " – TSS (edgeR ", nm, " sig overlay)"),
         pch = 1, cex = 0.4, col = rgb(0.5, 0.5, 0.5, 0.2),
         xlim = xlim_tss, ylim = ylim_tss)
    abline(h = 0, col = "black")
    abline(h = c(-1, 1), col = "grey60", lty = 3)

    if (any(both & fin_tss))
      points(r$tss$std_error[both & fin_tss], r$tss$estimate[both & fin_tss],
             pch = 19, cex = 0.7, col = "red")
    if (any(edger_only & fin_tss))
      points(r$tss$std_error[edger_only & fin_tss], r$tss$estimate[edger_only & fin_tss],
             pch = 19, cex = 0.7, col = "blue")

    legend("topleft",
           legend = c(paste0("Both (", sum(both), ")"),
                      paste0("edgeR only (", sum(edger_only), ")"),
                      paste0("Non-sig (", sum(!sig_edger), ")")),
           col = c("red", "blue", rgb(0.5, 0.5, 0.5, 0.4)),
           pch = c(19, 19, 1), cex = 0.7, bg = "white")

    # ---- RIGHT: edgeR MA plot — FILTERED coords, TSS sig overlay ----
    plot(filt_cpm[neither_k & fin_edg_f],
         filt_fc[neither_k & fin_edg_f],
         xlab = "logCPM (filtered)", ylab = "logFC (filtered)",
         main = paste0(comp, " – edgeR ", nm,
                       " FILTERED (", sum(keep), " OTUs, TSS sig overlay)"),
         pch = 1, cex = 0.4, col = rgb(0.5, 0.5, 0.5, 0.2),
         xlim = xlim_edg, ylim = ylim_edg)
    abline(h = 0, col = "black")
    abline(h = c(-1, 1), col = "grey60", lty = 3)

    if (any(both_k & fin_edg_f))
      points(filt_cpm[both_k & fin_edg_f],
             filt_fc[both_k & fin_edg_f],
             pch = 19, cex = 0.7, col = "red")
    if (any(tss_only_k & fin_edg_f))
      points(filt_cpm[tss_only_k & fin_edg_f],
             filt_fc[tss_only_k & fin_edg_f],
             pch = 19, cex = 0.7, col = "orange")

    legend("topleft",
           legend = c(paste0("Both (", sum(both_k), ")"),
                      paste0("TSS only (", sum(tss_only_k), ")"),
                      paste0("Non-sig (", sum(neither_k), ")")),
           col = c("red", "orange", rgb(0.5, 0.5, 0.5, 0.4)),
           pch = c(19, 19, 1), cex = 0.7, bg = "white")

    par(mfrow = c(1, 1))
  }
  dev.off()
  cat(sprintf("  Saved: TSS vs edgeR %s filtered effect cross plots\n", nm))
}

cat("\nDone. All filtered-coord plots saved to:\n  ", plot_dir, "\n")
