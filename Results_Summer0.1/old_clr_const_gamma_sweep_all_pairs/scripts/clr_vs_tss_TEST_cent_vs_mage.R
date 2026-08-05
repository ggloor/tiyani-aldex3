# =============================================================================
# TEST: CLR vs TSS (with constant) — cent vs mage only
# =============================================================================

library(ALDEx3)

# ---- Config -----------------------------------------------------------------
datadir   <- "~/Desktop/3383/0_git/tiyani-aldex3/data/"
gamma_val <- 1e-5
nsample   <- 128

# ---- Custom TSS scale with constant ----------------------------------------
tss_const <- function(X, logComp, gamma = 0.5, constant = 3) {
  P  <- nrow(X)
  ns <- dim(logComp)[3]
  LambdaScale <- matrix(rnorm(P * ns, constant, gamma), P, ns)
  logScale <- t(X) %*% LambdaScale
  return(logScale)
}

# ---- Load only cent and mage ------------------------------------------------
load(paste0(datadir, "cent.Rda"))
load(paste0(datadir, "mage.Rda"))

Y     <- cbind(cent, mage)
conds <- c(rep("cent", ncol(cent)), rep("mage", ncol(mage)))
data  <- data.frame(condition = conds)

# ---- Run ALDEx3 -------------------------------------------------------------
cat("Running CLR...\n")
res_clr <- ALDEx3::aldex(Y, ~condition, data,
                         nsample = nsample, scale = clr.sm, gamma = gamma_val)
sum_clr <- summary(res_clr)

cat("Running TSS...\n")
res_tss <- ALDEx3::aldex(Y, ~condition, data,
                         nsample = nsample, scale = tss_const, gamma = gamma_val)
sum_tss <- summary(res_tss)

# ---- Scatter plot -----------------------------------------------------------
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
     main = paste0("cent vs mage – CLR vs TSS(const=", formals(tss_const)$constant,
                   ") (gamma = ", gamma_val, ")"),
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

cat("Done! Both:", sum(both), " TSS only:", sum(tss_only),
    " CLR only:", sum(clr_only), " Neither:", sum(neither), "\n")
