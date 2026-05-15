# 03b_plots_other_norms.R
# Load saved .Rda results from 01b (unfiltered) and 02b (filtered),
# plot scatter comparisons for TMMwsp and RLE against ALDEx3 CLR and TSS.
# Inf values from -log10(0) are handled at plot time.

set.seed(12345)

# --- Paths and parameters ---
base_out  <- "~/Desktop/Western_COOP/Results_COOP/aldex_vs_edger"
gamma_val <- 1e-05
fdr_cut   <- 0.05

# --- Normalizations to plot (TMM already plotted by QuickStart) ---
edger_norms <- c("TMMwsp", "RLE")

# --- Cohort names and pairs ---
cohorts <- c("kin", "pup", "mid", "you", "mage", "eld", "cent")
pairs   <- combn(cohorts, 2)


# --- Scatter plot function (parameterized for any edgeR norm) ---
plot_scatter_matched <- function(aldex_padj, edger_fdr, aldex_label,
                                 edger_label, comp_name, n_total = NULL,
                                 filtered = FALSE) {
  
  x <- -log10(edger_fdr)
  y <- -log10(aldex_padj)
  
  finite <- is.finite(x) & is.finite(y)
  n_inf  <- sum(!finite)
  x <- x[finite]
  y <- y[finite]
  edger_fdr_f  <- edger_fdr[finite]
  aldex_padj_f <- aldex_padj[finite]
  
  subtitle <- sprintf("%d OTUs plotted", length(x))
  if (n_inf > 0) subtitle <- sprintf("%s (%d Inf removed)", subtitle, n_inf)
  if (!is.null(n_total)) subtitle <- sprintf("%s of %d", subtitle, n_total)
  filt_tag <- if (filtered) " filtered" else ""
  title <- sprintf("%s: %s vs %s%s (gamma = %g)\n%s",
                   comp_name, aldex_label, edger_label, filt_tag,
                   gamma_val, subtitle)
  
  sig_edger <- edger_fdr_f  < fdr_cut
  sig_aldex <- aldex_padj_f < fdr_cut
  
  n_both       <- sum(sig_edger & sig_aldex)
  n_edger_only <- sum(sig_edger & !sig_aldex)
  n_aldex_only <- sum(!sig_edger & sig_aldex)
  n_neither    <- sum(!sig_edger & !sig_aldex)
  
  cols <- rep("grey60", length(x))
  cols[sig_edger & sig_aldex]  <- "red"
  cols[sig_edger & !sig_aldex] <- "blue"
  cols[!sig_edger & sig_aldex] <- "orange"
  
  if (length(x) == 0) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "No OTUs to plot", cex = 1.2)
    return(invisible(NULL))
  }
  
  thr <- -log10(fdr_cut)
  lim <- max(c(x, y), na.rm = TRUE) * 1.1
  
  plot(x, y,
       xlab = bquote(-log[10](FDR) ~ .(edger_label) ~ "(edgeR)"),
       ylab = bquote(-log[10](FDR) ~ .(aldex_label) ~ "(ALDEx3)"),
       main = title,
       pch = 19, cex = 0.6, col = cols,
       xlim = c(0, lim), ylim = c(0, lim))
  
  abline(0, 1, col = "black", lty = 2)
  abline(h = thr, col = "grey40", lty = 3)
  abline(v = thr, col = "grey40", lty = 3)
  
  legend("topleft",
         legend = c(sprintf("Both (%d)", n_both),
                    sprintf("%s only (%d)", edger_label, n_edger_only),
                    sprintf("%s only (%d)", aldex_label, n_aldex_only),
                    sprintf("Neither (%d)", n_neither)),
         col = c("red", "blue", "orange", "grey60"),
         pch = 19, cex = 0.85, bty = "n")
  
  legend("bottomright",
         legend = c("FDR = 0.05", "y = x"),
         lty = c(3, 2), col = c("grey40", "black"),
         cex = 0.7, bty = "n", inset = c(0.02, 0.02))
}


# --- Load all 21 .Rda files from a directory ---
load_results <- function(dir_path) {
  results <- list()
  for (i in 1:ncol(pairs)) {
    fname <- paste0(pairs[1, i], "_vs_", pairs[2, i], ".Rda")
    load(file.path(dir_path, fname))
    results[[result$comparison]] <- result
  }
  return(results)
}


# =============================================================================
# UNFILTERED PLOTS
# =============================================================================

for (norm in edger_norms) {
  
  norm_dir <- file.path(base_out, paste0("unfiltered_", norm))
  dir.create(file.path(norm_dir, "plots"), showWarnings = FALSE)
  
  cat(sprintf("Loading unfiltered %s results...\n", norm))
  res <- load_results(norm_dir)
  
  # CLR vs this norm
  pdf(file.path(norm_dir, "plots",
                sprintf("clr_vs_%s_all_pairs.pdf", norm)), width = 7, height = 7)
  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- res[[comp]]
    plot_scatter_matched(r$clr$pval_adj, r$edger$FDR, "CLR", norm, comp)
  }
  dev.off()
  cat(sprintf("  Saved: CLR vs %s\n", norm))
  
  # TSS vs this norm
  pdf(file.path(norm_dir, "plots",
                sprintf("tss_vs_%s_all_pairs.pdf", norm)), width = 7, height = 7)
  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- res[[comp]]
    plot_scatter_matched(r$tss$pval_adj, r$edger$FDR, "TSS", norm, comp)
  }
  dev.off()
  cat(sprintf("  Saved: TSS vs %s\n", norm))
  
  # summary table
  cat(sprintf("\n=== Unfiltered %s: significant OTUs (FDR < 0.05) ===\n", norm))
  cat(sprintf("%-20s %6s %6s %8s %8s %8s\n",
              "Comparison", "CLR", "TSS", norm, "CLR+norm", "TSS+norm"))
  cat(paste(rep("-", 68), collapse = ""), "\n")
  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- res[[comp]]
    sig_clr  <- r$clr$pval_adj < fdr_cut
    sig_tss  <- r$tss$pval_adj < fdr_cut
    sig_norm <- r$edger$FDR    < fdr_cut
    cat(sprintf("%-20s %6d %6d %8d %8d %8d\n", comp,
                sum(sig_clr), sum(sig_tss), sum(sig_norm),
                sum(sig_clr & sig_norm), sum(sig_tss & sig_norm)))
  }
  cat("\n")
}


# =============================================================================
# FILTERED PLOTS
# =============================================================================

for (norm in edger_norms) {
  
  norm_dir <- file.path(base_out, paste0("filtered_", norm))
  dir.create(file.path(norm_dir, "plots"), showWarnings = FALSE)
  
  cat(sprintf("Loading filtered %s results...\n", norm))
  res <- load_results(norm_dir)
  
  # CLR vs this norm
  pdf(file.path(norm_dir, "plots",
                sprintf("clr_vs_%s_filtered_all_pairs.pdf", norm)), width = 7, height = 7)
  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- res[[comp]]
    plot_scatter_matched(r$clr$pval_adj, r$edger$FDR, "CLR", norm, comp,
                         n_total = r$n_kept, filtered = TRUE)
  }
  dev.off()
  cat(sprintf("  Saved: CLR vs %s filtered\n", norm))
  
  # TSS vs this norm
  pdf(file.path(norm_dir, "plots",
                sprintf("tss_vs_%s_filtered_all_pairs.pdf", norm)), width = 7, height = 7)
  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- res[[comp]]
    plot_scatter_matched(r$tss$pval_adj, r$edger$FDR, "TSS", norm, comp,
                         n_total = r$n_kept, filtered = TRUE)
  }
  dev.off()
  cat(sprintf("  Saved: TSS vs %s filtered\n", norm))
  
  # summary table
  cat(sprintf("\n=== Filtered %s: significant OTUs (FDR < 0.05) ===\n", norm))
  cat(sprintf("%-20s %8s %6s %6s %8s %8s %8s\n",
              "Comparison", "Kept", "CLR", "TSS", norm, "CLR+norm", "TSS+norm"))
  cat(paste(rep("-", 74), collapse = ""), "\n")
  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- res[[comp]]
    sig_clr  <- r$clr$pval_adj < fdr_cut
    sig_tss  <- r$tss$pval_adj < fdr_cut
    sig_norm <- r$edger$FDR    < fdr_cut
    cat(sprintf("%-20s %8d %6d %6d %8d %8d %8d\n", comp, r$n_kept,
                sum(sig_clr), sum(sig_tss), sum(sig_norm),
                sum(sig_clr & sig_norm), sum(sig_tss & sig_norm)))
  }
  cat("\n")
}

cat("All plots done.\n")