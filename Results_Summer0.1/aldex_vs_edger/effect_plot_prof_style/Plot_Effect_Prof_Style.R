# Plot_Effect_Prof_Style.R
# ============================================================================
# Single-panel ALDEx3 effect plots with layered significance overlay.
# Prof Gloor's preferred style:
#   - x-axis: std_error * sqrt(238)
#   - y-axis: estimate
#   - Base layer: all OTUs (black open circles)
#   - Layer 2: edgeR-significant OTUs (orange, pch=19)
#   - Layer 3: ALDEx-significant OTUs (red, pch=19, cex=0.5)
#
# If an OTU is significant in both, red (ALDEx) is drawn on top of orange
# (edgeR), so it appears red.
#
# Produces TWO versions per ALDEx norm × edgeR norm:
#   - Filtered:   edgeR significance from filtered run (filterByExpr + padded FDR)
#   - Unfiltered: edgeR significance from unfiltered run (all 1117 OTUs tested)
#
# 2 ALDEx norms × 3 edgeR norms × 2 filter versions = 12 PDFs, 21 pages each.
#
# Requires: .Rda files from Updated_Results0.1.R with fdr_unfilt field.
# ============================================================================


# #############################################################################
# SECTION 1: SETUP
# #############################################################################

rda_dir  <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/aldex_vs_edger/effect_cross_edgeR_allnorms"
out_dir  <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/aldex_vs_edger/effect_plot_prof_style"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cut      <- 0.05
gamma_val    <- 0.3
sqrt_n       <- sqrt(238)
norm_methods <- c("TMM", "TMMwsp", "RLE")

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

  load(fname)
  all_results[[comp]] <- result
  cat(sprintf("  Loaded: %s\n", comp))
}

cat(sprintf("\nLoaded %d comparisons.\n\n", length(all_results)))


# #############################################################################
# SECTION 3: PLOTTING FUNCTION
# #############################################################################
#
# One function handles all combos. Arguments:
#   aldex_norm  = "clr" or "tss"
#   edger_norm  = "TMM", "TMMwsp", or "RLE"
#   fdr_field   = "fdr" (filtered) or "fdr_unfilt" (unfiltered)
#   label       = string for title/filename ("filtered" or "unfiltered")
# #############################################################################

plot_prof_style <- function(aldex_norm, edger_norm, fdr_field, label) {

  aldex_label <- toupper(aldex_norm)  # "CLR" or "TSS"

  pdf_name <- sprintf("effect_%s_edgeR_%s_%s_gamma%s.pdf",
                       aldex_label, edger_norm, label,
                       gsub("\\.", "", as.character(gamma_val)))
  pdf(file.path(plot_dir, pdf_name), width = 8, height = 7)

  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- all_results[[comp]]

    # ALDEx data
    est <- r[[aldex_norm]]$estimate
    se  <- r[[aldex_norm]]$std_error
    x   <- se * sqrt_n
    y   <- est

    # Significance masks
    sig_ald  <- r[[aldex_norm]]$pval_adj    < fdr_cut
    sig_edge <- r$edger[[edger_norm]][[fdr_field]] < fdr_cut

    # Finite protection
    fin <- is.finite(x) & is.finite(y)

    # Axis limits with padding
    xlim <- range(x[fin], na.rm = TRUE) * c(0.95, 1.05)
    ylim <- range(y[fin], na.rm = TRUE) * 1.05

    # --- Base layer: all OTUs (black open circles) ---
    plot(x[fin], y[fin],
         xlab = paste0("Std Error × √238 (", aldex_label, ")"),
         ylab = paste0("Estimate (", aldex_label, ")"),
         main = paste0(comp, " – ", aldex_label,
                       " (edgeR ", edger_norm, " ", label, ")"),
         pch = 1, cex = 0.6,
         xlim = xlim, ylim = ylim)
    abline(h = 0, col = "black")

    # --- Layer 2: edgeR-significant (orange) ---
    if (any(sig_edge & fin))
      points(x[sig_edge & fin], y[sig_edge & fin],
             col = "orange", pch = 19)

    # --- Layer 3: ALDEx-significant (red, smaller, on top) ---
    if (any(sig_ald & fin))
      points(x[sig_ald & fin], y[sig_ald & fin],
             col = "red", pch = 19, cex = 0.5)

    # Legend
    n_both     <- sum(sig_ald & sig_edge)
    n_edge_only <- sum(sig_edge & !sig_ald)
    n_ald_only  <- sum(sig_ald & !sig_edge)
    n_neither   <- sum(!sig_ald & !sig_edge)

    legend("topleft",
           legend = c(paste0("Both (", n_both, ")"),
                      paste0("edgeR only (", n_edge_only, ")"),
                      paste0(aldex_label, " only (", n_ald_only, ")"),
                      paste0("Neither (", n_neither, ")")),
           col = c("red", "orange", "red", "black"),
           pch = c(19, 19, 19, 1),
           pt.cex = c(0.5, 1, 0.5, 0.6),
           cex = 0.8, bg = "white")
  }

  dev.off()
  cat(sprintf("  Saved: %s\n", pdf_name))
}


# #############################################################################
# SECTION 4: GENERATE ALL 12 PDFs
# #############################################################################

cat("===== Plotting — prof style =====\n\n")

for (nm in norm_methods) {
  cat(sprintf("--- edgeR %s ---\n", nm))

  # CLR × filtered
  plot_prof_style("clr", nm, "fdr",        "filtered")
  # CLR × unfiltered
  plot_prof_style("clr", nm, "fdr_unfilt", "unfiltered")

  # TSS × filtered
  plot_prof_style("tss", nm, "fdr",        "filtered")
  # TSS × unfiltered
  plot_prof_style("tss", nm, "fdr_unfilt", "unfiltered")
}

cat("\nDone. All plots saved to:\n  ", plot_dir, "\n")
