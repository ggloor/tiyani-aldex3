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
  eld  = eld,
  kin  = kin,
  mage = mage,
  mid  = mid,
  pup  = pup,
  you  = you
)

# Define ordered age groups (youngest to oldest)
group_order <- c("kin", "pup", "mid", "you", "mage", "eld", "cent")

# Group labels for plotting
group_labels <- c(
  kin  = "Kindergarten",
  pup  = "Primary",
  mid  = "Middle School",
  you  = "Youth",
  mage = "Middle Age",
  eld  = "Elderly",
  cent = "Centenarian"
)

# Get all pair combinations
pairs <- combn(names(datasets), 2)

# Create a data frame to store per-pair results
results <- data.frame(
  pair             = character(),
  mean_estimate    = numeric(),
  median_estimate  = numeric(),
  stringsAsFactors = FALSE
)

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
  data  <- data.frame(condition = conds)
  
  print(paste("Running:", name1, "vs", name2))
  
  # Run ALDEx3 with gamma = 0.3
  res_gamma <- ALDEx3::aldex(Y, ~condition, data, nsample=128, scale=tss.sm, gamma=0.3)
  sum_gamma <- summary(res_gamma)
  
  # Store mean and median estimates for both groups in this pair
  results <- rbind(results, data.frame(
    pair            = paste(name1, "vs", name2),
    group           = name1,
    mean_estimate   = mean(sum_gamma$estimate),
    median_estimate = median(sum_gamma$estimate),
    stringsAsFactors = FALSE
  ))
  results <- rbind(results, data.frame(
    pair            = paste(name1, "vs", name2),
    group           = name2,
    mean_estimate   = mean(sum_gamma$estimate),
    median_estimate = median(sum_gamma$estimate),
    stringsAsFactors = FALSE
  ))
}

# Summarise: average mean and median effect size per group across all comparisons
group_summary <- aggregate(cbind(mean_estimate, median_estimate) ~ group, data = results, FUN = mean)

# Order by age group
group_summary$group <- factor(group_summary$group, levels = group_order)
group_summary <- group_summary[order(group_summary$group), ]

print(group_summary)

# Save plots to PDF
pdf("~/Desktop/3383/0_git/tiyani-aldex3/Results/age_vs_effect_size_tss.pdf", width=8, height=6)

# ----- PLOT 1: Age Group vs Mean Effect Size -----
par(mar = c(7, 5, 4, 2))
plot(as.numeric(group_summary$group), group_summary$mean_estimate,
     main  = "Age Group vs Mean Effect Size (gamma = 0.3)",
     xlab  = "",
     ylab  = "Mean Effect Size",
     pch   = 19, col = "steelblue", cex = 1.4,
     xaxt  = "n",
     type  = "b",
     lty   = 2)
axis(1, at = 1:length(group_order),
     labels = group_labels[group_order],
     las = 2, cex.axis = 0.85)
abline(h = 0, lty = 2, col = "grey60")

# ----- PLOT 2: Age Group vs Median Effect Size -----
par(mar = c(7, 5, 4, 2))
plot(as.numeric(group_summary$group), group_summary$median_estimate,
     main  = "Age Group vs Median Effect Size (gamma = 0.3)",
     xlab  = "",
     ylab  = "Median Effect Size",
     pch   = 19, col = "darkorange", cex = 1.4,
     xaxt  = "n",
     type  = "b",
     lty   = 2)
axis(1, at = 1:length(group_order),
     labels = group_labels[group_order],
     las = 2, cex.axis = 0.85)
abline(h = 0, lty = 2, col = "grey60")

dev.off()