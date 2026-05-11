# 01_run_unfiltered.R
# Run ALDEx3 (CLR + TSS) and edgeR TMM on all 1117 OTUs for each pair.
# Saves one .Rda per comparison.

library(ALDEx3)
library(edgeR)
set.seed(12345)

# --- Paths and parameters ---
data_dir  <- "~/Desktop/Western_COOP/datacopy"
out_dir   <- "~/Desktop/Western_COOP/Results_COOP/aldex_vs_edger/unfiltered"
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
  
  # ALDEx3 CLR
  res_clr <- ALDEx3::aldex(Y, ~condition, dd,
                           nsample = nsample, scale = clr.sm, gamma = gamma_val)
  sum_clr <- summary(res_clr)
  
  # ALDEx3 TSS
  res_tss <- ALDEx3::aldex(Y, ~condition, dd,
                           nsample = nsample, scale = tss.sm, gamma = gamma_val)
  sum_tss <- summary(res_tss)
  
  # edgeR TMM (no filtering)
  group <- factor(conds)
  dge   <- DGEList(counts = Y, group = group)
  dge   <- calcNormFactors(dge, method = "TMM")
  dge   <- estimateDisp(dge)
  et    <- exactTest(dge)
  tt    <- topTags(et, n = Inf, sort.by = "none")$table
  
  # save
  result <- list(
    clr = data.frame(OTU = rownames(sum_clr),
                     estimate = sum_clr$estimate,
                     pval_adj = sum_clr$p.val.adj),
    tss = data.frame(OTU = rownames(sum_tss),
                     estimate = sum_tss$estimate,
                     pval_adj = sum_tss$p.val.adj),
    edger = data.frame(OTU = rownames(tt),
                       logFC = tt$logFC,
                       FDR   = tt$FDR),
    comparison = comp
  )
  
  fname <- paste0(name1, "_vs_", name2, ".Rda")
  save(result, file = file.path(out_dir, fname))
}

cat("Done. Results saved to:", out_dir, "\n")