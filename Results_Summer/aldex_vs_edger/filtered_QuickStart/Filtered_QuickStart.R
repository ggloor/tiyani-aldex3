# 02_run_filtered.R
# Same as 01 but edgeR uses filterByExpr first.
# ALDEx3 runs on the full OTU set (needs complete composition for normalization),
# then we subset its results to match the edgeR-kept OTUs.
# Inf rows from ALDEx3 MC sampling are removed before matching.

library(ALDEx3)
library(edgeR)
set.seed(12345)

# --- Paths and parameters ---
data_dir  <- "~/Desktop/Western_COOP/datacopy"
out_dir   <- "~/Desktop/Western_COOP/Results_COOP/aldex_vs_edger/filtered_QuickStart"
gamma_val <- 1e-05
nsample   <- 128

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

# --- Run all 21 comparisons ---
for (i in 1:ncol(pairs)) {
  
  name1 <- pairs[1, i]
  name2 <- pairs[2, i]
  comp  <- paste(name1, "vs", name2)
  cat(sprintf("[%2d/21] %s\n", i, comp))
  
  Y     <- cbind(datasets[[name1]], datasets[[name2]])
  conds <- c(rep(name1, ncol(datasets[[name1]])),
             rep(name2, ncol(datasets[[name2]])))
  dd    <- data.frame(condition = conds)
  
  # --- ALDEx3 CLR (full OTU set, needed for proper normalization) ---
  res_clr <- ALDEx3::aldex(Y, ~condition, dd,
                           nsample = nsample, scale = clr.sm, gamma = gamma_val)
  sum_clr <- summary(res_clr)
  
  # --- ALDEx3 TSS (full OTU set) ---
  res_tss <- ALDEx3::aldex(Y, ~condition, dd,
                           nsample = nsample, scale = tss.sm, gamma = gamma_val)
  sum_tss <- summary(res_tss)
  
  
  # --- edgeR QL pipeline with filterByExpr ---
  group  <- factor(conds)
  dge    <- DGEList(counts = Y, group = group)
  keep   <- filterByExpr(dge, group = group)
  cat(sprintf("  filterByExpr kept %d / %d OTUs\n", sum(keep), nrow(Y)))
  
  dge    <- dge[keep, , keep.lib.sizes = FALSE]
  dge    <- normLibSizes(dge)
  design <- model.matrix(~group)
  fit    <- glmQLFit(dge, design)
  qlf    <- glmQLFTest(fit, coef = 2)
  tt     <- topTags(qlf, n = Inf, sort.by = "none")$table
  
  # --- Subset ALDEx3 to match filterByExpr kept OTUs ---
  sum_clr_f <- sum_clr[keep, ]
  sum_tss_f <- sum_tss[keep, ]
  stopifnot(nrow(sum_clr_f) == nrow(tt))
  
  result <- list(
    clr = data.frame(OTU      = rownames(tt),
                     estimate = sum_clr_f$estimate,
                     pval_adj = sum_clr_f$p.val.adj),
    tss = data.frame(OTU      = rownames(tt),
                     estimate = sum_tss_f$estimate,
                     pval_adj = sum_tss_f$p.val.adj),
    edger = data.frame(OTU   = rownames(tt),
                       logFC = tt$logFC,
                       FDR   = tt$FDR),
    comparison = comp,
    n_kept = sum(keep)
  )
  
  fname <- paste0(name1, "_vs_", name2, ".Rda")
  save(result, file = file.path(out_dir, fname))
}

cat("Done. Results saved to:", out_dir, "\n")