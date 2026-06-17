# Plot_Pval_Scatter_AldexVsEdgeR.R
# ============================================================================
# -log10(FDR) scatter plots: ALDEx3 (CLR or TSS) vs edgeR (TMM/TMMwsp/RLE)
#
# Style: matches the CLR-vs-TSS reference plot (clr_vs_tss_(gmma=0.3)_plots.pdf)
#   - x-axis: -log10(FDR) edgeR
#   - y-axis: -log10(FDR) ALDEx3
#   - 1:1 diagonal (dashed), threshold lines at FDR=0.05 (dotted)
#   - Red = both sig, Blue = edgeR only, Orange = ALDEx only, Grey = neither
#
# Uses FILTERED OTUs only (those that passed filterByExpr).
# Zero p-values handled with min(nonzero)/10 substitution.
#
# Produces 6 PDFs (CLR/TSS x TMM/TMMwsp/RLE), each with 21 comparison pages.
#
# Requires: .Rda files from Updated_Results0.1.R (with coords_unfilt,
#           coords_filt, keep, fdr fields).
# ============================================================================


# #############################################################################
# SECTION 1: SETUP
# #############################################################################

rda_dir  <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/aldex_vs_edger/effect_cross_edgeR_allnorms"
out_dir  <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/aldex_vs_edger/pval_scatter_aldex_vs_edgeR(gamma0)"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cut      <- 0.05
gamma_val    <- 00
norm_methods <- c("TMM", "TMMwsp", "RLE")
aldex_norms  <- c("clr", "tss")

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
# SECTION 3: -log10(FDR) SCATTER PLOTS — ALDEx vs edgeR (FILTERED OTUs)
# #############################################################################
#
# For each ALDEx norm (CLR, TSS) x edgeR norm (TMM, TMMwsp, RLE):
#   - Subset to OTUs that passed filterByExpr (keep mask)
#   - x = -log10(edgeR FDR),  y = -log10(ALDEx adjusted p-value)
#   - Zero p-values replaced with min(nonzero)/10 before transform
#   - One page per comparison (21 pages per PDF)
# #############################################################################

cat("===== Generating -log10(FDR) scatter plots =====\n\n")

for (a_norm in aldex_norms) {
  for (e_norm in norm_methods) {

    a_label <- toupper(a_norm)   # "CLR" or "TSS"
    cat(sprintf("--- %s vs edgeR %s ---\n", a_label, e_norm))

    pdf_name <- sprintf("pval_scatter_%s_vs_edgeR_%s_filtered_gamma%s.pdf",
                        a_label, e_norm,
                        gsub("\\.", "", as.character(gamma_val)))
    pdf(file.path(plot_dir, pdf_name), width = 8, height = 7)

    for (i in 1:ncol(pairs)) {
      comp <- paste(pairs[1, i], "vs", pairs[2, i])
      r <- all_results[[comp]]

      # ---- Get the keep mask for this edgeR norm ----
      keep <- r$edger[[e_norm]]$keep

      # ---- Subset p-values to kept OTUs only ----
      aldex_padj <- r[[a_norm]]$pval_adj[keep]
      edger_fdr  <- r$edger[[e_norm]]$fdr[keep]
      n_kept     <- sum(keep)

      # ---- Zero p-value handling: min(nonzero)/10 substitution ----
      all_pvals   <- c(aldex_padj, edger_fdr)
      nonzero     <- all_pvals[all_pvals > 0]

      if (length(nonzero) > 0) {
        min_nonzero <- min(nonzero)
        aldex_padj[aldex_padj == 0] <- min_nonzero / 10
        edger_fdr[edger_fdr == 0]   <- min_nonzero / 10
      }

      # ---- -log10 transform ----
      x <- -log10(edger_fdr)
      y <- -log10(aldex_padj)
      threshold <- -log10(fdr_cut)

      # ---- Significance categories ----
      sig_edger <- x >= threshold
      sig_aldex <- y >= threshold

      n_both       <- sum(sig_edger & sig_aldex)
      n_edger_only <- sum(sig_edger & !sig_aldex)
      n_aldex_only <- sum(!sig_edger & sig_aldex)
      n_neither    <- sum(!sig_edger & !sig_aldex)

      # ---- Count zero p-values for subtitle ----
      n_zero <- sum(all_pvals == 0)

      # ---- Assign colors ----
      cols <- rep("grey60", length(x))
      cols[sig_edger & !sig_aldex] <- "blue"
      cols[!sig_edger & sig_aldex] <- "orange"
      cols[sig_edger & sig_aldex]  <- "red"

      # ---- Plot ----
      plot(x, y,
           xlab = paste0("-log10(FDR) edgeR ", e_norm),
           ylab = paste0("-log10(FDR) ", a_label),
           main = paste0(comp, " – ", a_label, " vs edgeR ", e_norm,
                         " (gamma = ", gamma_val, ")"),
           sub = if (n_zero > 0)
                   sprintf("%d filtered OTUs (%d zero p-vals replaced)",
                           n_kept, n_zero)
                 else
                   sprintf("%d filtered OTUs", n_kept),
           pch = 19, cex = 0.6, col = cols)

      abline(0, 1, col = "black", lty = 2)         # 1:1 diagonal
      abline(h = threshold, col = "grey40", lty = 3) # horizontal FDR threshold
      abline(v = threshold, col = "grey40", lty = 3) # vertical FDR threshold

      legend("topleft",
             legend = c(
               paste0("Both (", n_both, ")"),
               paste0("edgeR only (", n_edger_only, ")"),
               paste0(a_label, " only (", n_aldex_only, ")"),
               paste0("Neither (", n_neither, ")")
             ),
             col = c("red", "blue", "orange", "grey60"),
             pch = 19, cex = 0.8)
    }

    dev.off()
    cat(sprintf("  Saved: %s\n", pdf_name))
  }
}


# #############################################################################
# SECTION 4: SUMMARY TABLE — sig counts per comparison (filtered OTUs)
# #############################################################################

cat("\n===== Summary: significance overlap (filtered OTUs) =====\n\n")

for (a_norm in aldex_norms) {
  for (e_norm in norm_methods) {

    a_label <- toupper(a_norm)
    cat(sprintf("=== %s vs edgeR %s (FDR < %.2f, gamma = %g, filtered) ===\n\n",
                a_label, e_norm, fdr_cut, gamma_val))
    cat(sprintf("%-18s %6s %6s %6s %8s %10s %10s\n",
                "Comparison", "Kept", a_label, "edgeR",
                "Both", paste0(a_label, " only"), "edgeR only"))
    cat(paste(rep("-", 78), collapse = ""), "\n")

    for (i in 1:ncol(pairs)) {
      comp <- paste(pairs[1, i], "vs", pairs[2, i])
      r <- all_results[[comp]]

      keep       <- r$edger[[e_norm]]$keep
      aldex_padj <- r[[a_norm]]$pval_adj[keep]
      edger_fdr  <- r$edger[[e_norm]]$fdr[keep]

      sig_aldex <- aldex_padj < fdr_cut
      sig_edger <- edger_fdr  < fdr_cut

      cat(sprintf("%-18s %6d %6d %6d %8d %10d %10d\n",
                  comp, sum(keep),
                  sum(sig_aldex), sum(sig_edger),
                  sum(sig_aldex & sig_edger),
                  sum(sig_aldex & !sig_edger),
                  sum(sig_edger & !sig_aldex)))
    }
    cat("\n")
  }
}

# #############################################################################
# SECTION 5: SAVE SUMMARY DATA TO CSV
# #############################################################################

cat("===== Saving summary CSV =====\n\n")

summary_rows <- list()
row_idx <- 0

for (a_norm in aldex_norms) {
  for (e_norm in norm_methods) {
    for (i in 1:ncol(pairs)) {
      comp <- paste(pairs[1, i], "vs", pairs[2, i])
      r <- all_results[[comp]]

      keep       <- r$edger[[e_norm]]$keep
      aldex_padj <- r[[a_norm]]$pval_adj[keep]
      edger_fdr  <- r$edger[[e_norm]]$fdr[keep]

      sig_aldex <- aldex_padj < fdr_cut
      sig_edger <- edger_fdr  < fdr_cut

      row_idx <- row_idx + 1
      summary_rows[[row_idx]] <- data.frame(
        aldex_norm  = toupper(a_norm),
        edger_norm  = e_norm,
        comparison  = comp,
        n_kept      = sum(keep),
        n_sig_aldex = sum(sig_aldex),
        n_sig_edger = sum(sig_edger),
        n_both      = sum(sig_aldex & sig_edger),
        n_aldex_only = sum(sig_aldex & !sig_edger),
        n_edger_only = sum(sig_edger & !sig_aldex),
        stringsAsFactors = FALSE
      )
    }
  }
}

summary_df <- do.call(rbind, summary_rows)
csv_path   <- file.path(out_dir, "pval_scatter_summary.csv")
write.csv(summary_df, csv_path, row.names = FALSE)
cat(sprintf("  Saved: %s\n", csv_path))

cat("\nDone. All plots and data saved to:\n  ", out_dir, "\n")
