
library(ALDEx3)

# Load all datasets
load('~/Desktop/3383/0_git/tiyani-aldex3/data/cent.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/eld.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/kin.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mage.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mid.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/pup.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/you.Rda')

# Store datasets in a list
datasets <- list(
  cent = cent,
  eld = eld,
  kin = kin,
  mage = mage,
  mid = mid,
  pup = pup,
  you = you
)

# Get all pair combinations
pairs <- combn(names(datasets), 2)
all_results <- list()
# START saving to PDF
pdf("~/Desktop/3383/0_git/tiyani-aldex3/Results/pairwise_analysis_plots_with_clr.pdf", width=8, height=6)
num <-0
# Loop through each pair
for (i in 1:ncol(pairs)) {
  
  name1 <- pairs[1, i]
  name2 <- pairs[2, i]
  
  data1 <- datasets[[name1]]
  data2 <- datasets[[name2]]
  
  # Combine data
  Y <- cbind(data1, data2)
  
  # Create condition labels
  conds <- c(rep(name1, ncol(data1)), rep(name2, ncol(data2)))
  data <- data.frame(condition = conds)
  
  print(paste("Running:", name1, "vs", name2))
  num <-print(num+1)
  
  # Run WITHOUT gamma (gamma ≈ 0)
  res_no_gamma <- ALDEx3::aldex(Y, ~condition, data, nsample=128, scale=clr.sm, gamma=1e-3)
  sum_no_gamma <- summary(res_no_gamma)
  
  # Run WITH gamma (gamma = 0.3)
  res_gamma <- ALDEx3::aldex(Y, ~condition, data, nsample=128, scale=clr.sm, gamma=0.3)
  sum_gamma <- summary(res_gamma)
  
  # Setting Sig values according to no gamma result
  sig <- which(sum_no_gamma$p.val.adj < 0.05)

  
  # Plot: With gamma
  plot(sum_gamma$std.error, sum_gamma$estimate,
       main = paste(name1, "vs", name2),
       sub = "With Gamma (gamma = 0.3)",
       xlab = "Std Error", ylab = "Estimate", col="darkgreen")
  sig_gamma <- which(sum_gamma$p.val.adj < 0.05)
  points(sum_gamma$std.error[sig], sum_gamma$estimate[sig], 
         pch=19, col='orange', cex=0.5)
  points(sum_gamma$std.error[sig_gamma], sum_gamma$estimate[sig_gamma], 
         pch=19, col='red', cex=0.5)
  abline(h=0)
}

all_results[[paste(name1, "vs", name2)]] <- data.frame(
  taxon = rownames(sum_gamma),
  estimate_gamma = sum_gamma$estimate,
  std_error_gamma = sum_gamma$std.error,
  pval_gamma = sum_gamma$p.val,
  pval_adj_gamma = sum_gamma$p.val.adj,
  estimate_no_gamma = sum_no_gamma$estimate,
  pval_adj_no_gamma = sum_no_gamma$p.val.adj,
  comparison = paste(name1, "vs", name2)
)


dev.off()

print("All plots saved to PDF!")