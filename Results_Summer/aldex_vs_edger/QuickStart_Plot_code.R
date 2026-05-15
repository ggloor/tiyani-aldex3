# 03_plots.R
# Load saved .Rda results from 01 and 02, make scatter plots.
# Inf values from -log10(0) are handled at plot time (not removed from data).

set.seed(12345)

# --- Paths and parameters ---
base_out   <- "~/Desktop/Western_COOP/Results_COOP/aldex_vs_edger"
unfilt_dir <- file.path(base_out, "unfiltered_QuickStart")
filt_dir   <- file.path(base_out, "filtered_QuickStart")
gamma_val  <- 1e-05
fdr_cut    <- 0.05

# --- Cohort names and pairs ---
cohorts <- c("kin", "pup", "mid", "you", "mage", "eld", "cent")
pairs   <- combn(cohorts, 2)

dir.create(file.path(unfilt_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(filt_dir,   "plots"), showWarnings = FALSE)


# --- Scatter plot function (matched significant OTUs only) ---
# Restricts to OTUs significant in BOTH ALDEx3 and edgeR (FDR < 0.05).
# Shows how the two methods compare on effect magnitude for shared hits.
plot_scatter_matched <- function(aldex_padj, edger_fdr, label, comp_name,
                                 n_total = NULL) {
  
  # -log10 transform
  x <- -log10(edger_fdr)
  y <- -log10(aldex_padj)
  
  # remove Inf values (where raw p-value was exactly 0)
  finite <- is.finite(x) & is.finite(y)
  n_inf  <- sum(!finite)
  x <- x[finite]
  y <- y[finite]
  edger_fdr_f  <- edger_fdr[finite]
  aldex_padj_f <- aldex_padj[finite]
  
  # build title with accurate count of plotted points
  subtitle <- sprintf("%d OTUs plotted", length(x))
  if (n_inf > 0) subtitle <- sprintf("%s (%d Inf removed)", subtitle, n_inf)
  if (!is.null(n_total)) subtitle <- sprintf("%s of %d", subtitle, n_total)
  title <- sprintf("%s: %s vs TMM (gamma = %g)\n%s", comp_name, label,
                   gamma_val, subtitle)
  
  # classify significance on the finite set
  sig_tmm <- edger_fdr_f  < fdr_cut
  sig_alx <- aldex_padj_f < fdr_cut
  
  # counts for legend (computed before subsetting)
  n_both     <- sum(sig_tmm & sig_alx)
  n_tmm_only <- sum(sig_tmm & !sig_alx)
  n_alx_only <- sum(!sig_tmm & sig_alx)
  
  # color by category
  cols <- rep("grey60", length(x))
  cols[sig_tmm & sig_alx]  <- "red"
  cols[sig_tmm & !sig_alx] <- "blue"
  cols[!sig_tmm & sig_alx] <- "orange"
  
  n_neither  <- sum(!sig_tmm & !sig_alx)
  
  if (length(x) == 0) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "No significant OTUs in either method", cex = 1.2)
    return(invisible(NULL))
  }
  
  thr <- -log10(fdr_cut)
  lim <- max(c(x, y), na.rm = TRUE) * 1.1
  
  plot(x, y,
       xlab = expression(-log[10](FDR) ~ "TMM (edgeR)"),
       ylab = bquote(-log[10](FDR) ~ .(label) ~ "(ALDEx3)"),
       main = title,
       pch = 19, cex = 0.6, col = cols,
       xlim = c(0, lim), ylim = c(0, lim))
  
  abline(0, 1, col = "black", lty = 2)
  abline(h = thr, col = "grey40", lty = 3)
  abline(v = thr, col = "grey40", lty = 3)
  
  legend("topleft",
         legend = c(sprintf("Both (%d)", n_both),
                    sprintf("TMM only (%d)", n_tmm_only),
                    sprintf("%s only (%d)", label, n_alx_only),
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

cat("Loading unfiltered results...\n")
res_uf <- load_results(unfilt_dir)

# CLR vs TMM
pdf(file.path(unfilt_dir, "plots", "clr_vs_tmm_QuickStart_all_pairs0.3.pdf"), width = 7, height = 7)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_uf[[comp]]
  plot_scatter_matched(r$clr$pval_adj, r$edger$FDR, "CLR", comp)
}
dev.off()
cat("Saved: unfiltered CLR vs TMM plots\n")

# TSS vs TMM
pdf(file.path(unfilt_dir, "plots", "tss_vs_tmm_QuickStart_all_pairs0.3.pdf"), width = 7, height = 7)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_uf[[comp]]
  plot_scatter_matched(r$tss$pval_adj, r$edger$FDR, "TSS", comp)
}
dev.off()
cat("Saved: unfiltered TSS vs TMM plots\n")

# summary
cat("\n=== Unfiltered: significant OTUs (FDR < 0.05) ===\n")
cat(sprintf("%-20s %6s %6s %6s %8s %8s\n",
            "Comparison", "CLR", "TSS", "TMM", "CLR+TMM", "TSS+TMM"))
cat(paste(rep("-", 64), collapse = ""), "\n")
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_uf[[comp]]
  sig_clr <- r$clr$pval_adj < fdr_cut
  sig_tss <- r$tss$pval_adj < fdr_cut
  sig_tmm <- r$edger$FDR    < fdr_cut
  cat(sprintf("%-20s %6d %6d %6d %8d %8d\n", comp,
              sum(sig_clr), sum(sig_tss), sum(sig_tmm),
              sum(sig_clr & sig_tmm), sum(sig_tss & sig_tmm)))
}


# =============================================================================
# FILTERED PLOTS
# =============================================================================

cat("\nLoading filtered results...\n")
res_f <- load_results(filt_dir)

# CLR vs TMM
pdf(file.path(filt_dir, "plots", "clr_vs_tmm_QuickStart_filtered_all_pairs0.3.pdf"), width = 7, height = 7)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_f[[comp]]
  plot_scatter_matched(r$clr$pval_adj, r$edger$FDR, "CLR", comp,
                       n_total = r$n_kept)
}
dev.off()
cat("Saved: filtered CLR vs TMM plots\n")

# TSS vs TMM
pdf(file.path(filt_dir, "plots", "tss_vs_tmm_QuickStart_filtered_all_pairs0.3.pdf"), width = 7, height = 7)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_f[[comp]]
  plot_scatter_matched(r$tss$pval_adj, r$edger$FDR, "TSS", comp,
                       n_total = r$n_kept)
}
dev.off()
cat("Saved: filtered TSS vs TMM plots\n")

# summary
cat("\n=== Filtered: significant OTUs (FDR < 0.05) ===\n")
cat(sprintf("%-20s %8s %6s %6s %6s %8s %8s\n",
            "Comparison", "Matched", "CLR", "TSS", "TMM", "CLR+TMM", "TSS+TMM"))
cat(paste(rep("-", 72), collapse = ""), "\n")
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_f[[comp]]
  sig_clr <- r$clr$pval_adj < fdr_cut
  sig_tss <- r$tss$pval_adj < fdr_cut
  sig_tmm <- r$edger$FDR    < fdr_cut
  cat(sprintf("%-20s %8d %6d %6d %6d %8d %8d\n", comp, r$n_kept,
              sum(sig_clr), sum(sig_tss), sum(sig_tmm),
              sum(sig_clr & sig_tmm), sum(sig_tss & sig_tmm)))
}

cat("\nAll plots done.\n")