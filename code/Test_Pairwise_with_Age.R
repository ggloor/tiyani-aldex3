
# Single Pairwise ALDEx3 Analysis - Test Run

library(ALDEx2)
library(ALDEx3)


# Load datasets (count data)

load('~/Desktop/3383/0_git/tiyani-aldex3/data/kin.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/pup.Rda')


# Load age metadata

load('~/Desktop/3383/0_git/tiyani-aldex3/data/kin.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/pup.age.Rda')


# Set up comparison: kin vs pup


# Combine count data
Y <- cbind(kin, pup)

# make a conditions vector
# for a continuous variable (age) we replace this with the ages of the samples
age <- c(as.numeric(kin.age[1,]), as.numeric(pup.age[1,]))
conds <- c(rep("K", length(kin.age)), rep('P', length(pup.age)))
data <- data.frame(age = age, condition = conds)

# Check counts
print("Condition counts:")
print(table(conds))

# Set seed for reproducibility
set.seed(42)


# Run ALDEx3 with gamma = 0.3

print("Running ALDEx3 analysis...")

res_gamma <- ALDEx3::aldex(Y, ~condition, data, nsample = 128, scale = clr.sm, gamma = 0.3)
sum_gamma <- summary(res_gamma)

# Get significant taxa
sig_gamma <- which(sum_gamma$p.val.adj < 0.05)

print(paste("Significant taxa:", length(sig_gamma)))

# Plot

plot(sum_gamma$std.error, sum_gamma$estimate,
     main = "K vs P (Kindergarten vs Primary)",
     sub = paste0("Gamma = 0.3 | n_sig = ", length(sig_gamma)),
     xlab = "Std Error", ylab = "Estimate", 
     col = "darkgreen", pch = 1)

# Highlight significant points
points(sum_gamma$std.error[sig_gamma], sum_gamma$estimate[sig_gamma], 
       pch = 19, col = 'red', cex = 0.5)
abline(h = 0, lty = 2)

print("Done!")