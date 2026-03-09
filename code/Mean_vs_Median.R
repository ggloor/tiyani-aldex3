library(ALDEx2)
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

# Define midpoint ages for each group
ages <- c(
  cent = 97,    # Centenarians (94+)
  eld = 70,     # Elderly (60-79)
  mage = 40,    # Middle-aged (30-50)
  you = 22,     # Youth (19-24)
  mid = 13.5,   # Middle school (13-14)
  pup = 10,     # Primary school (8-12)
  kin = 4.5     # Kindergarten (3-6)
)

# Get all pair combinations
pairs <- combn(names(datasets), 2)

# Create a data frame to store results
results <- data.frame(
  pair = character(),
  age_diff = numeric(),
  median_estimate = numeric(),
  mean_estimate = numeric(),
  mean_median_diff = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each pair
for (i in 1:ncol(pairs)) {
  
  name1 <- pairs[1, i]
  name2 <- pairs[2, i]
  
  data1 <- datasets[[name1]]
  data2 <- datasets[[name2]]
  
  # Calculate age difference (absolute value)
  age_diff <- abs(ages[name1] - ages[name2])
  
  # Combine data
  Y <- cbind(data1, data2)
  
  # Create condition labels
  conds <- c(rep(name1, ncol(data1)), rep(name2, ncol(data2)))
  data <- data.frame(condition = conds)
  
  print(paste("Running:", name1, "vs", name2, "| Age diff:", age_diff))
  
  # Run WITH gamma (gamma = 0.3)
  res_gamma <- ALDEx3::aldex(Y, ~condition, data, nsample=128, scale=clr.sm, gamma=0.3)
  sum_gamma <- summary(res_gamma)
  
  # Calculate median, mean, and difference
  med <- median(sum_gamma$estimate)
  mn <- mean(sum_gamma$estimate)
  diff <- mn - med
  
  # Store results
  results <- rbind(results, data.frame(
    pair = paste(name1, "vs", name2),
    age_diff = age_diff,
    median_estimate = med,
    mean_estimate = mn,
    mean_median_diff = diff
  ))
}

# View results table
print(results)

# ----- PLOT -----

par(mar = c(5, 5, 4, 2))

plot(results$age_diff, results$mean_median_diff,
     main = "Mean-Median Difference vs Age Difference",
     xlab = "Age Difference (years)",
     ylab = "Mean - Median Effect Size",
     pch = 19, col = "blue", cex = 1.2,
     xlim = c(-5, 100),
     ylim = c(-0.5, 0.4))
abline(h = 0, lty = 2)
abline(lm(mean_median_diff ~ age_diff, data = results), col = "red", lwd = 2)

# Add labels
text(results$age_diff, results$mean_median_diff, 
     labels = results$pair, 
     pos = 4,
     cex = 0.55,
     offset = 0.3)

