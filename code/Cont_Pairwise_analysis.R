
library(ALDEx3)

# Load count data
load('~/Desktop/3383/0_git/tiyani-aldex3/data/cent.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/eld.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/kin.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mage.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mid.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/pup.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/you.Rda')

# Load age data
load('~/Desktop/3383/0_git/tiyani-aldex3/data/cent.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/eld.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/kin.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mage.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mid.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/pup.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/you.age.Rda')

# Combine all count data
Y_all <- cbind(kin, pup, mid, you, mage, eld, cent)

# Combine all ages
all_ages <- c(as.numeric(kin.age[1,]), as.numeric(pup.age[1,]),
              as.numeric(mid.age[1,]), as.numeric(you.age[1,]),
              as.numeric(mage.age[1,]), as.numeric(eld.age[1,]),
              as.numeric(cent.age[1,]))

data_all <- data.frame(age = all_ages)

# Run continuous model across full age range
print("Running continuous model across all age groups...")
res <- ALDEx3::aldex(Y_all, ~age, data_all, nsample=128, scale=clr.sm, gamma=0.3)
s <- summary(res)

# Results
results <- data.frame(
  taxon = s$entity,
  estimate = s$estimate,
  std_error = s$std.error,
  pval_adj = s$p.val.adj
)

write.csv(results, "~/Desktop/3383/0_git/tiyani-aldex3/Results/continuous_all_groups_clr.csv", row.names = FALSE)

# Volcano plot
pdf("~/Desktop/3383/0_git/tiyani-aldex3/Results/continuous_all_groups_plots.pdf", width=8, height=6)

cols <- ifelse(s$p.val.adj < 0.05 & s$estimate > 0, "red",
               ifelse(s$p.val.adj < 0.05 & s$estimate < 0, "blue", "grey60"))

plot(s$estimate, -log10(s$p.val.adj),
     xlab = "Estimate (CLR change per year of age)",
     ylab = "-log10(adjusted p-value)",
     main = "Continuous Age Model - All Groups (gamma = 0.3)",
     pch = 19, cex = 0.6, col = cols)
abline(h = -log10(0.05), col = "grey40", lty = 3)
abline(v = 0, col = "grey40", lty = 3)

legend("topright",
       legend = c(
         paste0("Increase with age (", sum(cols == "red"), ")"),
         paste0("Decrease with age (", sum(cols == "blue"), ")"),
         paste0("Not significant (", sum(cols == "grey60"), ")")
       ),
       col = c("red", "blue", "grey60"),
       pch = 19, cex = 0.7)

# Effect plot
cols2 <- ifelse(s$p.val.adj < 0.05, "red", "darkgreen")

plot(s$std.error, s$estimate,
     xlab = "Standard Error",
     ylab = "Estimate (CLR change per year of age)",
     main = "Effect Plot - Continuous Age Model (gamma = 0.3)",
     pch = 19, cex = 0.6, col = cols2)
abline(h = 0, col = "grey40")

legend("topright",
       legend = c(
         paste0("Significant (", sum(cols2 == "red"), ")"),
         paste0("Not significant (", sum(cols2 == "darkgreen"), ")")
       ),
       col = c("red", "darkgreen"),
       pch = 19, cex = 0.7)

dev.off()

print(paste("Total taxa:", nrow(s)))
print(paste("Significant taxa (FDR < 0.05):", sum(s$p.val.adj < 0.05)))
