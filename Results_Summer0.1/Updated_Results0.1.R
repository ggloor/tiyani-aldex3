# effect_cross_aldex_edgeR_allnorms.R
# Effect plots comparing ALDEx3 (CLR/TSS) vs edgeR (TMM, TMMwsp, RLE) with
# cross-method significance overlay.
#
# ALDEx3 runs on all 1117 OTUs (unfiltered).
# edgeR stores BOTH:
#   - coords_unfilt: logFC/logCPM from unfiltered run (all OTUs, for plotting all points)
#   - coords_filt:   logFC/logCPM from filtered run (kept OTUs only, different norm factors)
#   - fdr:           FDR from filtered run (non-kept OTUs padded with 1)
#   - fdr_unfilt:    FDR from unfiltered run (all 1117 OTUs tested)
# Gamma = 0.3.

library(ALDEx3)
library(edgeR)
set.seed(12345)

# --- Paths and parameters ---
data_dir <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/data"
out_dir  <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/aldex_vs_edger/effect_cross_edgeR_allnorms"
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

nsample   <- 128
gamma_val <- 0.3
fdr_cut   <- 0.05

# edgeR normalization methods to compare
norm_methods <- c("TMM", "TMMwsp", "RLE")

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
  design <- model.matrix(~group)
  
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
  
  # --- edgeR: loop over normalization methods ---
  edger_list <- list()
  
  for (nm in norm_methods) {
    cat(sprintf("  edgeR %s unfiltered (coordinates)\n", nm))
    dge_all <- DGEList(counts = Y, group = group)
    dge_all <- normLibSizes(dge_all, method = nm)
    fit_all <- glmQLFit(dge_all, design)
    qlf_all <- glmQLFTest(fit_all, coef = 2)
    tt_all  <- topTags(qlf_all, n = Inf, sort.by = "none")$table
    
    cat(sprintf("  edgeR %s filtered (FDR)\n", nm))
    dge_f <- DGEList(counts = Y, group = group)
    keep  <- filterByExpr(dge_f, group = group)
    cat(sprintf("  filterByExpr kept %d / %d OTUs\n", sum(keep), nrow(Y)))
    
    dge_f <- dge_f[keep, , keep.lib.sizes = FALSE]
    dge_f <- normLibSizes(dge_f, method = nm)
    fit_f <- glmQLFit(dge_f, design)
    qlf_f <- glmQLFTest(fit_f, coef = 2)
    tt_f  <- topTags(qlf_f, n = Inf, sort.by = "none")$table
    
    edger_fdr_full <- rep(1, nrow(Y))
    names(edger_fdr_full) <- rownames(Y)
    edger_fdr_full[keep] <- tt_f$FDR
    
    # ---- Store BOTH unfiltered and filtered coordinates ----
    # coords_unfilt: logFC/logCPM/F/SE from the unfiltered run (all OTUs)
    #   - used for plotting all OTUs on MA plot or effect-space plot
    # coords_filt:   logFC/logCPM/F/SE from the filtered run (kept OTUs only)
    #   - different norm factors (only kept OTUs contribute)
    # F:             QL F-statistic from glmQLFTest
    # SE:            standard error of logFC, derived as |logFC| / sqrt(F)
    # keep:          logical mask — which OTUs passed filterByExpr
    # fdr:           FDR from filtered run, non-kept OTUs padded with 1
    # fdr_unfilt:    FDR from unfiltered run (all 1117 OTUs tested, no filterByExpr)
    edger_list[[nm]] <- list(
      coords_unfilt = data.frame(OTU    = rownames(tt_all),
                                 logFC  = tt_all$logFC,
                                 logCPM = tt_all$logCPM,
                                 F_stat = tt_all$F,
                                 SE     = ifelse(tt_all$F > 0,
                                                 abs(tt_all$logFC) / sqrt(tt_all$F),
                                                 NA)),
      coords_filt   = data.frame(OTU    = rownames(tt_f),
                                 logFC  = tt_f$logFC,
                                 logCPM = tt_f$logCPM,
                                 F_stat = tt_f$F,
                                 SE     = ifelse(tt_f$F > 0,
                                                 abs(tt_f$logFC) / sqrt(tt_f$F),
                                                 NA)),
      keep       = keep,
      fdr        = edger_fdr_full,
      fdr_unfilt = tt_all$FDR,
      n_kept     = sum(keep)
    )
  }
  
  # --- Save ---
  result <- list(
    clr = data.frame(OTU       = rownames(sum_clr),
                     estimate  = sum_clr$estimate,
                     std_error = sum_clr$std.error,
                     pval_adj  = sum_clr$p.val.adj),
    tss = data.frame(OTU       = rownames(sum_tss),
                     estimate  = sum_tss$estimate,
                     std_error = sum_tss$std.error,
                     pval_adj  = sum_tss$p.val.adj),
    edger      = edger_list,
    comparison = comp
  )
  all_results[[comp]] <- result
  
  fname <- paste0(name1, "_vs_", name2, ".Rda")
  save(result, file = file.path(out_dir, fname))
}

cat("\nAll results saved to:", out_dir, "\n\n")


# #############################################################################
# PART 2: PLOTTING \u2014 Prof Gloor style
# #############################################################################
#
# Single-panel ALDEx3 effect plots with layered significance overlay.
#   x-axis: std_error * sqrt(238)
#   y-axis: estimate
#   Colour scheme (all 4 categories distinct):
#     - Neither:    grey open circles (base layer)
#     - edgeR only: orange filled (pch=19)
#     - ALDEx only: blue filled (pch=19, cex=0.5)
#     - Both:       red filled (pch=19, cex=0.5, drawn last / on top)
#
# Produces BOTH filtered and unfiltered versions in separate folders:
#   plots_filtered/   \u2014 edgeR FDR from filterByExpr run
#   plots_unfiltered/ \u2014 edgeR FDR from all-OTU run (no filterByExpr)
#
# 2 ALDEx norms \u00d7 3 edgeR norms \u00d7 2 filter versions = 12 PDFs.
# #############################################################################

cat("===== PART 2: PLOTS (prof style) =====\n\n")

sqrt_n <- sqrt(238)

# Create separate output folders for filtered and unfiltered
plot_dir_filt   <- file.path(out_dir, "plots_filtered")
plot_dir_unfilt <- file.path(out_dir, "plots_unfiltered")
dir.create(plot_dir_filt,   recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir_unfilt, recursive = TRUE, showWarnings = FALSE)

# --- Helper function: one PDF of 21 pages ---
plot_prof_style <- function(aldex_norm, edger_norm, fdr_field, label,
                            plot_out_dir) {

  aldex_label <- toupper(aldex_norm)

  pdf_name <- sprintf("effect_%s_edgeR_%s_%s_gamma%s.pdf",
                       aldex_label, edger_norm, label,
                       gsub("\\.", "", as.character(gamma_val)))
  pdf(file.path(plot_out_dir, pdf_name), width = 8, height = 7)

  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- all_results[[comp]]

    # ALDEx data
    est <- r[[aldex_norm]]$estimate
    se  <- r[[aldex_norm]]$std_error
    x   <- se * sqrt_n
    y   <- est

    # Significance masks
    sig_ald  <- r[[aldex_norm]]$pval_adj          < fdr_cut
    sig_edge <- r$edger[[edger_norm]][[fdr_field]] < fdr_cut

    both       <- sig_ald & sig_edge
    edge_only  <- sig_edge & !sig_ald
    ald_only   <- sig_ald & !sig_edge
    neither    <- !sig_ald & !sig_edge

    # Finite protection
    fin <- is.finite(x) & is.finite(y)

    # Axis limits
    xlim <- range(x[fin], na.rm = TRUE) * c(0.95, 1.05)
    ylim <- range(y[fin], na.rm = TRUE) * 1.05

    # Layer 1: neither \u2014 grey open circles (base)
    plot(x[neither & fin], y[neither & fin],
         xlab = expression("Std Error" %*% sqrt(238)),
         ylab = "Estimate",
         main = paste0(comp, " \u2013 ", aldex_label,
                       " (edgeR ", edger_norm, " ", label, ")"),
         pch = 1, cex = 0.6, col = "grey60",
         xlim = xlim, ylim = ylim)
    abline(h = 0, col = "black")

    # Layer 2: edgeR only \u2014 orange
    if (any(edge_only & fin))
      points(x[edge_only & fin], y[edge_only & fin],
             col = "orange", pch = 19)

    # Layer 3: ALDEx only \u2014 blue
    if (any(ald_only & fin))
      points(x[ald_only & fin], y[ald_only & fin],
             col = "blue", pch = 19, cex = 0.5)

    # Layer 4: Both \u2014 red (on top of everything)
    if (any(both & fin))
      points(x[both & fin], y[both & fin],
             col = "red", pch = 19, cex = 0.5)

    # Legend
    legend("topleft",
           legend = c(paste0("Both (", sum(both), ")"),
                      paste0("edgeR only (", sum(edge_only), ")"),
                      paste0(aldex_label, " only (", sum(ald_only), ")"),
                      paste0("Neither (", sum(neither), ")")),
           col = c("red", "orange", "blue", "grey60"),
           pch = c(19, 19, 19, 1),
           pt.cex = c(0.5, 1, 0.5, 0.6),
           cex = 0.8, bg = "white")
  }

  dev.off()
  cat(sprintf("  Saved: %s \u2192 %s\n", pdf_name, label))
}

# --- Generate all 12 PDFs ---
for (nm in norm_methods) {
  cat(sprintf("--- edgeR %s ---\n", nm))
  plot_prof_style("clr", nm, "fdr",        "filtered",   plot_dir_filt)
  plot_prof_style("clr", nm, "fdr_unfilt", "unfiltered", plot_dir_unfilt)
  plot_prof_style("tss", nm, "fdr",        "filtered",   plot_dir_filt)
  plot_prof_style("tss", nm, "fdr_unfilt", "unfiltered", plot_dir_unfilt)
}


# #############################################################################
# PART 3: SUMMARY TABLE
# #############################################################################

cat("\n===== PART 3: SUMMARY =====\n\n")

for (nm in norm_methods) {
  
  cat(sprintf("=== ALDEx3 vs edgeR %s (FDR < 0.05, gamma = %g) ===\n\n",
              nm, gamma_val))
  cat(sprintf("%-18s %6s %6s %6s %8s %8s %9s %9s\n",
              "Comparison", "Kept", "CLR", "TSS", "edgeR",
              "CLR+edgR", "TSS+edgR", "edgR only"))
  cat(paste(rep("-", 85), collapse = ""), "\n")
  
  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- all_results[[comp]]
    
    sig_clr   <- r$clr$pval_adj    < fdr_cut
    sig_tss   <- r$tss$pval_adj    < fdr_cut
    sig_edger <- r$edger[[nm]]$fdr < fdr_cut
    
    cat(sprintf("%-18s %6d %6d %6d %6d %8d %8d %9d\n", comp,
                r$edger[[nm]]$n_kept,
                sum(sig_clr), sum(sig_tss), sum(sig_edger),
                sum(sig_clr & sig_edger),
                sum(sig_tss & sig_edger),
                sum(sig_edger & !sig_clr)))
  }
  cat("\n")
}

cat("Done.\n")