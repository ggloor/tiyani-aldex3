# effect_cross_aldex_edgeR_allnorms.R
# Effect plots comparing ALDEx3 (CLR/TSS) vs edgeR (TMM, TMMwsp, RLE) with
# cross-method significance overlay.
#
# ALDEx3 runs on all 1117 OTUs (unfiltered).
# edgeR runs unfiltered for coordinates, filtered for FDR.
# Gamma = 0.3.

library(ALDEx3)
library(edgeR)
set.seed(12345)

# --- Paths and parameters ---
data_dir <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/data"
out_dir  <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer/aldex_vs_edger/effect_cross_edgeR_allnorms"
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
    
    edger_list[[nm]] <- list(
      coords = data.frame(OTU    = rownames(tt_all),
                          logFC  = tt_all$logFC,
                          logCPM = tt_all$logCPM),
      fdr    = edger_fdr_full,
      n_kept = sum(keep)
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
# PART 2: PLOTTING
# #############################################################################

cat("===== PART 2: PLOTS =====\n\n")

# Loop over each normalization method, produce CLR vs edgeR and TSS vs edgeR PDFs
for (nm in norm_methods) {
  
  # =========================================================================
  # PLOT SET: CLR vs edgeR (current norm)
  # =========================================================================
  cat(sprintf("--- CLR vs edgeR %s ---\n", nm))
  
  pdf(file.path(out_dir, "plots",
                sprintf("effect_cross_CLR_edgeR_%s_gamma%s.pdf",
                        nm, gsub("\\.", "", as.character(gamma_val)))),
      width = 12, height = 6)
  
  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- all_results[[comp]]
    
    sig_clr   <- r$clr$pval_adj    < fdr_cut
    sig_edger <- r$edger[[nm]]$fdr < fdr_cut
    
    both       <- sig_clr & sig_edger
    edger_only <- sig_edger & !sig_clr
    clr_only   <- sig_clr & !sig_edger
    neither    <- !sig_clr & !sig_edger
    
    xlim_clr <- range(r$clr$std_error, na.rm = TRUE) * c(0.95, 1.05)
    xlim_edg <- range(r$edger[[nm]]$coords$logCPM, na.rm = TRUE) * c(0.95, 1.05)
    ylim_clr <- range(r$clr$estimate,  na.rm = TRUE) * 1.05
    ylim_edg <- range(r$edger[[nm]]$coords$logFC, na.rm = TRUE) * 1.05
    
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    
    # ---- Left: CLR effect plot, edgeR sig overlay ----
    plot(r$clr$std_error[neither], r$clr$estimate[neither],
         xlab = "Std Error", ylab = "Estimate",
         main = paste0(comp, " \u2013 CLR (edgeR ", nm, " sig overlay)"),
         pch = 1, cex = 0.4, col = rgb(0.5, 0.5, 0.5, 0.2),
         xlim = xlim_clr, ylim = ylim_clr)
    abline(h = 0, col = "black")
    abline(h = c(-1, 1), col = "grey60", lty = 3)
    
    if (any(both))
      points(r$clr$std_error[both], r$clr$estimate[both],
             pch = 19, cex = 0.7, col = "red")
    if (any(edger_only))
      points(r$clr$std_error[edger_only], r$clr$estimate[edger_only],
             pch = 19, cex = 0.7, col = "blue")
    
    legend("topleft",
           legend = c(paste0("Both (", sum(both), ")"),
                      paste0("edgeR only (", sum(edger_only), ")"),
                      paste0("Non-sig (", sum(!sig_edger), ")")),
           col = c("red", "blue", rgb(0.5, 0.5, 0.5, 0.4)),
           pch = c(19, 19, 1), cex = 0.7, bg = "white")
    
    # ---- Right: edgeR MA plot, CLR sig overlay ----
    plot(r$edger[[nm]]$coords$logCPM[neither],
         r$edger[[nm]]$coords$logFC[neither],
         xlab = "logCPM", ylab = "logFC",
         main = paste0(comp, " \u2013 edgeR ", nm, " (CLR sig overlay)"),
         pch = 1, cex = 0.4, col = rgb(0.5, 0.5, 0.5, 0.2),
         xlim = xlim_edg, ylim = ylim_edg)
    abline(h = 0, col = "black")
    abline(h = c(-1, 1), col = "grey60", lty = 3)
    
    if (any(both))
      points(r$edger[[nm]]$coords$logCPM[both],
             r$edger[[nm]]$coords$logFC[both],
             pch = 19, cex = 0.7, col = "red")
    if (any(clr_only))
      points(r$edger[[nm]]$coords$logCPM[clr_only],
             r$edger[[nm]]$coords$logFC[clr_only],
             pch = 19, cex = 0.7, col = "orange")
    
    legend("topleft",
           legend = c(paste0("Both (", sum(both), ")"),
                      paste0("CLR only (", sum(clr_only), ")"),
                      paste0("Non-sig (", sum(!sig_clr), ")")),
           col = c("red", "orange", rgb(0.5, 0.5, 0.5, 0.4)),
           pch = c(19, 19, 1), cex = 0.7, bg = "white")
    
    par(mfrow = c(1, 1))
  }
  dev.off()
  cat(sprintf("  Saved: CLR vs edgeR %s effect cross plots\n", nm))
  
  
  # =========================================================================
  # PLOT SET: TSS vs edgeR (current norm)
  # =========================================================================
  cat(sprintf("--- TSS vs edgeR %s ---\n", nm))
  
  pdf(file.path(out_dir, "plots",
                sprintf("effect_cross_TSS_edgeR_%s_gamma%s.pdf",
                        nm, gsub("\\.", "", as.character(gamma_val)))),
      width = 12, height = 6)
  
  for (i in 1:ncol(pairs)) {
    comp <- paste(pairs[1, i], "vs", pairs[2, i])
    r <- all_results[[comp]]
    
    sig_tss   <- r$tss$pval_adj    < fdr_cut
    sig_edger <- r$edger[[nm]]$fdr < fdr_cut
    
    both       <- sig_tss & sig_edger
    edger_only <- sig_edger & !sig_tss
    tss_only   <- sig_tss & !sig_edger
    neither    <- !sig_tss & !sig_edger
    
    xlim_tss <- range(r$tss$std_error, na.rm = TRUE) * c(0.95, 1.05)
    xlim_edg <- range(r$edger[[nm]]$coords$logCPM, na.rm = TRUE) * c(0.95, 1.05)
    ylim_tss <- range(r$tss$estimate,  na.rm = TRUE) * 1.05
    ylim_edg <- range(r$edger[[nm]]$coords$logFC, na.rm = TRUE) * 1.05
    
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    
    # ---- Left: TSS effect plot, edgeR sig overlay ----
    plot(r$tss$std_error[neither], r$tss$estimate[neither],
         xlab = "Std Error", ylab = "Estimate",
         main = paste0(comp, " \u2013 TSS (edgeR ", nm, " sig overlay)"),
         pch = 1, cex = 0.4, col = rgb(0.5, 0.5, 0.5, 0.2),
         xlim = xlim_tss, ylim = ylim_tss)
    abline(h = 0, col = "black")
    abline(h = c(-1, 1), col = "grey60", lty = 3)
    
    if (any(both))
      points(r$tss$std_error[both], r$tss$estimate[both],
             pch = 19, cex = 0.7, col = "red")
    if (any(edger_only))
      points(r$tss$std_error[edger_only], r$tss$estimate[edger_only],
             pch = 19, cex = 0.7, col = "blue")
    
    legend("topleft",
           legend = c(paste0("Both (", sum(both), ")"),
                      paste0("edgeR only (", sum(edger_only), ")"),
                      paste0("Non-sig (", sum(!sig_edger), ")")),
           col = c("red", "blue", rgb(0.5, 0.5, 0.5, 0.4)),
           pch = c(19, 19, 1), cex = 0.7, bg = "white")
    
    # ---- Right: edgeR MA plot, TSS sig overlay ----
    plot(r$edger[[nm]]$coords$logCPM[neither],
         r$edger[[nm]]$coords$logFC[neither],
         xlab = "logCPM", ylab = "logFC",
         main = paste0(comp, " \u2013 edgeR ", nm, " (TSS sig overlay)"),
         pch = 1, cex = 0.4, col = rgb(0.5, 0.5, 0.5, 0.2),
         xlim = xlim_edg, ylim = ylim_edg)
    abline(h = 0, col = "black")
    abline(h = c(-1, 1), col = "grey60", lty = 3)
    
    if (any(both))
      points(r$edger[[nm]]$coords$logCPM[both],
             r$edger[[nm]]$coords$logFC[both],
             pch = 19, cex = 0.7, col = "red")
    if (any(tss_only))
      points(r$edger[[nm]]$coords$logCPM[tss_only],
             r$edger[[nm]]$coords$logFC[tss_only],
             pch = 19, cex = 0.7, col = "orange")
    
    legend("topleft",
           legend = c(paste0("Both (", sum(both), ")"),
                      paste0("TSS only (", sum(tss_only), ")"),
                      paste0("Non-sig (", sum(!sig_tss), ")")),
           col = c("red", "orange", rgb(0.5, 0.5, 0.5, 0.4)),
           pch = c(19, 19, 1), cex = 0.7, bg = "white")
    
    par(mfrow = c(1, 1))
  }
  dev.off()
  cat(sprintf("  Saved: TSS vs edgeR %s effect cross plots\n", nm))
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