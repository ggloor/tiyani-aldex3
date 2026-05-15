library(ALDEx2)
library(ALDEx3)

# load Count data
load('~/Desktop/3383/0_git/tiyani-aldex3/data/cent.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/eld.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/kin.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mage.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mid.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/pup.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/you.Rda')



# load Age metadata
load('~/Desktop/3383/0_git/tiyani-aldex3/data/cent.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/eld.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/kin.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mage.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mid.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/pup.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/you.age.Rda')



# Organize Data into lists
counts_list <- list(
  kin=kin, 
  pup=pup,
  mid=mid, 
  you=you,
  mage=mage,
  eld=eld,
  cent=cent
)

age_list <- list(
  kin=kin.age, 
  pup=pup.age,
  mid=mid.age, 
  you=you.age, 
  mage=mage.age, 
  eld=eld.age,
  cent=cent.age
)



# Validate the Data
cat("Data Validation \n")
for (name in names(counts_list)) { 
  n_counts <- ncol(counts_list[[name]])
  n_age <- ncol(age_list[[name]])
  status <- ifelse(n_counts == n_age, "ok", "MisMatch") 
  cat(name, "-counts:", n_counts, " , Age:", n_age, status,"\n")
  }
cat("\n")



# Setup & Pairwise Analyses

group_names <- names(counts_list)
pairs <- combn(group_names, 2 , simplify = FALSE)

#Start saving to PDf
#pdf("")

for (i in seq_along(pairs)){ 
  pair <- pairs[[i]]
  g1 <- pair[1]
  g2 <- pair[2]
  pair_name <- paste( g1, "vs", g2, seq = " ")

  cat("Running Comparison", i, "of", length (pairs), ":", pair_name, "\n")
  
  Y <- cbind(counts_list[[g1]], counts_list[[g2]])
  
  age <- c(as.numeric(age_list[[g1]][1,]), as.numeric(age_list[[g2]][1,]))
  conds <- c(rep(g1, length(age_list[[g1]])), rep(g2, length(age_list[[g2]])))
  data <- data.frame(age = age, condition = conds)
  
  # Run Analyses
  res_age <- ALDEx3::aldex(Y, ~age, data, nsmaple=128, scale=clr.sm, gamma=0.3)
  sum_age <- summary(res_age)
  
  sig_gamma <- which(sum_gamma$p.val.adj < 0.05)
  
  print(paste("Significant taxa:", length(sig_gamma)))
  
 
  
  # Plot 1: Condition (Discrete)
  plot(sum_conds$std.error, sum_conds$estimate,
       main = paste(title_name, "- Condition"),
       sub = paste0("Gamma = 0.3 | n_sig = ", length(sig_conds)),
       xlab = "Std Error", ylab = "Estimate",
       col = "darkgreen", pch = 1)
  points(sum_conds$std.error[sig_conds], sum_conds$estimate[sig_conds],
         pch = 19, col = "red", cex = 0.5)
  abline(h = 0, lty = 2)
  legend("topright", 
         legend = c("Non-significant", "Significant (p.adj < 0.05)"),
         col = c("darkgreen", "red"),
         pch = c(1, 19),
         bty = "n")
  
  # Plot 2: Age (Continuous)
  plot(sum_age$std.error, sum_age$estimate,
       main = paste(title_name, "- Age"),
       sub = paste0("Gamma = 0.3 | n_sig = ", length(sig_age)),
       xlab = "Std Error", ylab = "Estimate",
       col = "darkgreen", pch = 1)
  points(sum_age$std.error[sig_age], sum_age$estimate[sig_age],
         pch = 19, col = "red", cex = 0.5)
  abline(h = 0, lty = 2)
  legend("topright", 
         legend = c("Non-significant", "Significant (p.adj < 0.05)"),
         col = c("darkgreen", "red"),
         pch = c(1, 19),
         bty = "n")
  
  }

#dev.off()
