# 03_plots.R
# Load saved .Rda results from 01 and 02, make all scatter plots.
# FDR-only classification (no effect size), matching the CLR vs TSS style.

set.seed(12345)

# --- Paths and parameters ---
base_out   <- "~/Desktop/Western_COOP/Results_COOP/aldex_vs_edger"
unfilt_dir <- file.path(base_out, "unfiltered")
filt_dir   <- file.path(base_out, "filtered")
gamma_val  <- 1e-05
fdr_cut    <- 0.05

# --- Cohort names and pairs (same order as 01/02) ---
cohorts <- c("kin", "pup", "mid", "you", "mage", "eld", "cent")
pairs   <- combn(cohorts, 2)

dir.create(file.path(unfilt_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(filt_dir,   "plots"), showWarnings = FALSE)


# --- Scatter plot function ---
# just FDR-based coloring, plain filled dots, like the CLR vs TSS plots
plot_scatter <- function(aldex_padj, edger_fdr, label, title) {
  
  x <- -log10(edger_fdr)
  y <- -log10(aldex_padj)
  thr <- -log10(fdr_cut)
  
  # cap infinities
  cap <- max(c(x[is.finite(x)], y[is.finite(y)]), na.rm = TRUE) * 1.1
  x[!is.finite(x)] <- cap
  y[!is.finite(y)] <- cap
  
  # color by FDR only
  sig_tmm <- edger_fdr < fdr_cut
  sig_alx <- aldex_padj < fdr_cut
  
  cols <- rep("grey60", length(x))
  cols[sig_tmm & !sig_alx] <- "blue"
  cols[!sig_tmm & sig_alx] <- "orange"
  cols[sig_tmm & sig_alx]  <- "red"
  
  # counts
  n_both <- sum(sig_tmm & sig_alx)
  n_tmm  <- sum(sig_tmm & !sig_alx)
  n_alx  <- sum(!sig_tmm & sig_alx)
  n_none <- sum(!sig_tmm & !sig_alx)
  
  plot(x, y,
       xlab = expression(-log[10](FDR) ~ "TMM (edgeR)"),
       ylab = bquote(-log[10](FDR) ~ .(label) ~ "(ALDEx3)"),
       main = title,
       pch = 19, cex = 0.6, col = cols,
       xlim = c(0, cap), ylim = c(0, cap))
  
  abline(0, 1, col = "black", lty = 2)
  abline(h = thr, col = "grey40", lty = 3)
  abline(v = thr, col = "grey40", lty = 3)
  
  legend("topleft",
         legend = c(sprintf("Both (%d)", n_both),
                    sprintf("TMM only (%d)", n_tmm),
                    sprintf("%s only (%d)", label, n_alx),
                    sprintf("Neither (%d)", n_none)),
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
pdf(file.path(unfilt_dir, "plots", "clr_vs_tmm_all_pairs.pdf"), width = 7, height = 7)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_uf[[comp]]
  plot_scatter(r$clr$pval_adj, r$edger$FDR,
               "CLR", paste0(comp, ": CLR vs TMM (gamma = ", gamma_val, ")"))
}
dev.off()
cat("Saved: unfiltered CLR vs TMM plots\n")

# TSS vs TMM
pdf(file.path(unfilt_dir, "plots", "tss_vs_tmm_all_pairs.pdf"), width = 7, height = 7)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_uf[[comp]]
  plot_scatter(r$tss$pval_adj, r$edger$FDR,
               "TSS", paste0(comp, ": TSS vs TMM (gamma = ", gamma_val, ")"))
}
dev.off()
cat("Saved: unfiltered TSS vs TMM plots\n")

# summary
cat("\n=== Unfiltered: significant OTUs (FDR < 0.05) ===\n")
cat(sprintf("%-20s %6s %6s %6s\n", "Comparison", "CLR", "TSS", "TMM"))
cat(paste(rep("-", 42), collapse = ""), "\n")
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_uf[[comp]]
  cat(sprintf("%-20s %6d %6d %6d\n", comp,
              sum(r$clr$pval_adj < fdr_cut),
              sum(r$tss$pval_adj < fdr_cut),
              sum(r$edger$FDR < fdr_cut)))
}


# =============================================================================
# FILTERED PLOTS
# =============================================================================

cat("\nLoading filtered results...\n")
res_f <- load_results(filt_dir)

# CLR vs TMM
pdf(file.path(filt_dir, "plots", "clr_vs_tmm_filtered_all_pairs.pdf"), width = 7, height = 7)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_f[[comp]]
  plot_scatter(r$clr$pval_adj, r$edger$FDR,
               "CLR", paste0(comp, ": CLR vs TMM filtered (gamma = ", gamma_val, ")\n",
                             r$n_kept, " OTUs retained by filterByExpr"))
}
dev.off()
cat("Saved: filtered CLR vs TMM plots\n")

# TSS vs TMM
pdf(file.path(filt_dir, "plots", "tss_vs_tmm_filtered_all_pairs.pdf"), width = 7, height = 7)
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_f[[comp]]
  plot_scatter(r$tss$pval_adj, r$edger$FDR,
               "TSS", paste0(comp, ": TSS vs TMM filtered (gamma = ", gamma_val, ")\n",
                             r$n_kept, " OTUs retained by filterByExpr"))
}
dev.off()
cat("Saved: filtered TSS vs TMM plots\n")

# summary
cat("\n=== Filtered: significant OTUs (FDR < 0.05) ===\n")
cat(sprintf("%-20s %8s %6s %6s %6s\n", "Comparison", "Retained", "CLR", "TSS", "TMM"))
cat(paste(rep("-", 52), collapse = ""), "\n")
for (i in 1:ncol(pairs)) {
  comp <- paste(pairs[1, i], "vs", pairs[2, i])
  r <- res_f[[comp]]
  cat(sprintf("%-20s %8d %6d %6d %6d\n", comp, r$n_kept,
              sum(r$clr$pval_adj < fdr_cut),
              sum(r$tss$pval_adj < fdr_cut),
              sum(r$edger$FDR < fdr_cut)))
}

cat("\nAll plots done.\n")