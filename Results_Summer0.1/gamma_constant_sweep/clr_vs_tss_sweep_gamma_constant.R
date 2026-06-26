# =============================================================================
# Sweep: CLR vs TSS across gamma and constant combinations
# All 21 pairwise comparisons x 5 gammas x 4 constants = 420 plots
# =============================================================================

library(ALDEx3)

# ---- Config -----------------------------------------------------------------
datadir <- "~/Desktop/3383/0_git/tiyani-aldex3/data/"
outdir  <- "~/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/"
nsample <- 128

gammas    <- c(1e-5, 1, 2, 3)
constants <- c(4)

# ---- Create output folder structure ----------------------------------------
sweepdir <- paste0(outdir, "gamma_constant_sweep/")
dir.create(sweepdir, showWarnings = FALSE, recursive = TRUE)
for (g in gammas) {
  dir.create(paste0(sweepdir, "gamma_", g, "/"), showWarnings = FALSE)
}

# ---- Load datasets ----------------------------------------------------------
load(paste0(datadir, "cent.Rda"))
load(paste0(datadir, "eld.Rda"))
load(paste0(datadir, "kin.Rda"))
load(paste0(datadir, "mage.Rda"))
load(paste0(datadir, "mid.Rda"))
load(paste0(datadir, "pup.Rda"))
load(paste0(datadir, "you.Rda"))

datasets <- list(
  cent = cent, eld = eld, kin = kin,
  mage = mage, mid = mid, pup = pup, you = you
)

pairs <- combn(names(datasets), 2)

# ---- Sweep ------------------------------------------------------------------
for (g in gammas) {
  for (k in constants) {

    k_label <- gsub("\\.", "p", as.character(k))
    pdfname <- paste0(sweepdir, "gamma_", g, "/clr_vs_tss_const", k_label, ".pdf")

    cat("=== gamma =", g, " | constant =", k, " ===\n")

    # Define custom TSS inside the loop so it captures current k
    tss_const <- function(X, logComp, gamma = 0.5) {
      P  <- nrow(X)
      ns <- dim(logComp)[3]
      LambdaScale <- matrix(rnorm(P * ns, k, gamma), P, ns)
      logScale <- t(X) %*% LambdaScale
      return(logScale)
    }

    pdf(pdfname, width = 8, height = 7)

    for (i in 1:ncol(pairs)) {

      name1 <- pairs[1, i]
      name2 <- pairs[2, i]
      comp_name <- paste(name1, "vs", name2)
      cat("  ", comp_name, "\n")

      Y     <- cbind(datasets[[name1]], datasets[[name2]])
      conds <- c(rep(name1, ncol(datasets[[name1]])),
                 rep(name2, ncol(datasets[[name2]])))
      data  <- data.frame(condition = conds)

      # CLR (standard)
      res_clr <- ALDEx3::aldex(Y, ~condition, data,
                               nsample = nsample, scale = clr.sm, gamma = g)
      sum_clr <- summary(res_clr)

      # TSS (with constant)
      res_tss <- ALDEx3::aldex(Y, ~condition, data,
                               nsample = nsample, scale = tss_const, gamma = g)
      sum_tss <- summary(res_tss)

      # p-values
      tss_padj <- sum_tss$p.val.adj
      clr_padj <- sum_clr$p.val.adj

      # Replace zero p-values with min(nonzero)/10
      all_pvals <- c(tss_padj, clr_padj)
      min_nz    <- min(all_pvals[all_pvals > 0])
      tss_padj[tss_padj == 0] <- min_nz / 10
      clr_padj[clr_padj == 0] <- min_nz / 10

      x <- -log10(tss_padj)
      y <- -log10(clr_padj)
      threshold <- -log10(0.05)

      sig_tss <- x >= threshold
      sig_clr <- y >= threshold

      both     <- sig_tss & sig_clr
      tss_only <- sig_tss & !sig_clr
      clr_only <- !sig_tss & sig_clr
      neither  <- !sig_tss & !sig_clr

      cols <- rep("grey60", length(x))
      cols[tss_only] <- "blue"
      cols[clr_only] <- "orange"
      cols[both]     <- "red"

      plot(x, y,
           xlab = expression(-log[10](FDR) ~ TSS),
           ylab = expression(-log[10](FDR) ~ CLR),
           main = paste0(comp_name, " – CLR vs TSS(const=", k,
                         ") (gamma = ", g, ")"),
           pch  = 19, cex = 0.6, col = cols)

      abline(0, 1, col = "black", lty = 2)
      abline(h = threshold, col = "grey40", lty = 3)
      abline(v = threshold, col = "grey40", lty = 3)

      legend("topleft",
             legend = c(paste0("Both (",     sum(both),     ")"),
                        paste0("TSS only (", sum(tss_only), ")"),
                        paste0("CLR only (", sum(clr_only), ")"),
                        paste0("Neither (",  sum(neither),  ")")),
             col = c("red", "blue", "orange", "grey60"),
             pch = 19, cex = 0.8)
    }

    dev.off()
    cat("  Saved:", pdfname, "\n\n")
  }
}

cat("All done! 20 PDFs saved.\n")
