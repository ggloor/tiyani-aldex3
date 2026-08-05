# =============================================================================
# Effect Cross Plots: CLR (with constant) vs TSS
# Sweep across gamma and constant combinations
# =============================================================================

library(ALDEx3)

# ---- Config -----------------------------------------------------------------
datadir <- "~/Desktop/3383/0_git/tiyani-aldex3/data/"
outdir  <- "/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results_Summer0.1/gamma_constant_sweep/Constant_effect_space"
nsample <- 128

gammas    <- c(1e-3)
constants <- seq(-1, 1, by = 0.1)

# ---- Create output folder structure ----------------------------------------
sweepdir <- paste0(outdir, "effect_cross_clr_constant_sweep/")
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
    pdfname <- paste0(sweepdir, "gamma_", g, "/effect_cross_clr_const", k_label, ".pdf")

    cat("=== gamma =", g, " | constant =", k, " ===\n")

    # Define custom CLR inside the loop so it captures current k
    clr_const <- function(X, logComp, gamma = 0.5) {
      P  <- nrow(X)
      ns <- dim(logComp)[3]
      LambdaScale <- matrix(rnorm(P * ns, k, gamma), P, ns)
      logScale <- t(X) %*% LambdaScale
      return(logScale)
    }

    pdf(pdfname, width = 14, height = 7)

    for (i in 1:ncol(pairs)) {

      name1 <- pairs[1, i]
      name2 <- pairs[2, i]
      comp_name <- paste(name1, "vs", name2)
      cat("  ", comp_name, "\n")

      Y     <- cbind(datasets[[name1]], datasets[[name2]])
      conds <- c(rep(name1, ncol(datasets[[name1]])),
                 rep(name2, ncol(datasets[[name2]])))
      data  <- data.frame(condition = conds)

      # CLR (with constant)
      res_clr <- ALDEx3::aldex(Y, ~condition, data,
                               nsample = nsample, scale = clr_const, gamma = g)
      sum_clr <- summary(res_clr)

      # TSS (standard)
      res_tss <- ALDEx3::aldex(Y, ~condition, data,
                               nsample = nsample, scale = tss.sm, gamma = g)
      sum_tss <- summary(res_tss)

      # Significance
      sig_clr <- sum_clr$p.val.adj < 0.05
      sig_tss <- sum_tss$p.val.adj < 0.05

      both     <- sig_clr & sig_tss
      tss_only <- sig_tss & !sig_clr
      clr_only <- sig_clr & !sig_tss
      neither  <- !sig_clr & !sig_tss

      # ---- Side-by-side panels ------------------------------------------------
      par(mfrow = c(1, 2))

      # --- Left panel: CLR effect space, TSS sig overlay ---
      cols_left <- rep("grey80", nrow(sum_clr))
      cols_left[both]     <- "red"
      cols_left[tss_only] <- "blue"

      plot(sum_clr$std.error, sum_clr$estimate,
           xlab = "Std Error", ylab = "Estimate",
           main = paste0(comp_name, " – CLR(const=", k, ") (TSS sig overlay)"),
           pch = 1, cex = 0.5, col = cols_left,
           font.main = 2)

      # Overlay significant points as filled
      if (any(both)) {
        points(sum_clr$std.error[both], sum_clr$estimate[both],
               col = "red", pch = 19, cex = 0.5)
      }
      if (any(tss_only)) {
        points(sum_clr$std.error[tss_only], sum_clr$estimate[tss_only],
               col = "blue", pch = 19, cex = 0.5)
      }

      abline(h = 0)
      abline(h = c(-1, 1), col = "grey60", lty = 3)

      legend("topleft",
             legend = c(paste0("Both (", sum(both), ")"),
                        paste0("TSS only (", sum(tss_only), ")"),
                        paste0("Non-sig (", sum(neither), ")")),
             col = c("red", "blue", "grey80"),
             pch = 19, cex = 0.8)

      # --- Right panel: TSS effect space, CLR sig overlay ---
      cols_right <- rep("grey80", nrow(sum_tss))
      cols_right[both]     <- "red"
      cols_right[clr_only] <- "orange"

      plot(sum_tss$std.error, sum_tss$estimate,
           xlab = "Std Error", ylab = "Estimate",
           main = paste0(comp_name, " – TSS (CLR(const=", k, ") sig overlay)"),
           pch = 1, cex = 0.5, col = cols_right,
           font.main = 2)

      # Overlay significant points as filled
      if (any(both)) {
        points(sum_tss$std.error[both], sum_tss$estimate[both],
               col = "red", pch = 19, cex = 0.5)
      }
      if (any(clr_only)) {
        points(sum_tss$std.error[clr_only], sum_tss$estimate[clr_only],
               col = "orange", pch = 19, cex = 0.5)
      }

      abline(h = 0)
      abline(h = c(-1, 1), col = "grey60", lty = 3)

      legend("topleft",
             legend = c(paste0("Both (", sum(both), ")"),
                        paste0("CLR only (", sum(clr_only), ")"),
                        paste0("Non-sig (", sum(neither), ")")),
             col = c("red", "orange", "grey80"),
             pch = 19, cex = 0.8)

      par(mfrow = c(1, 1))
    }

    dev.off()
    cat("  Saved:", pdfname, "\n\n")
  }
}

cat("All done!\n")
