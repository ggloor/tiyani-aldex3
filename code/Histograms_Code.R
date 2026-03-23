pvals <- read.csv("/Users/pranavdivvela/Desktop/3383/0_git/tiyani-aldex3/Results/pairwise_pvalues_clr.csv")

# ---- A: One histogram per comparison ----
pdf("~/Desktop/3383/0_git/tiyani-aldex3/Results/effect_size_histograms.pdf", width=8, height=6)

for (comp in unique(pvals$comparison)) {
  sub <- pvals[pvals$comparison == comp, ]
  
  hist(sub$estimate, breaks = 50,
       main = comp,
       xlab = "Estimate (effect size)",
       col = "steelblue", border = "white")
  abline(v = 0, col = "red", lty = 2, lwd = 2)
  abline(v = median(sub$estimate), col = "orange", lty = 2, lwd = 2)
  abline(v = mean(sub$estimate), col = "darkgreen", lty = 2, lwd = 2)
  
  legend("topright",
         legend = c(
           "Zero",
           paste0("Median (", round(median(sub$estimate), 3), ")"),
           paste0("Mean (", round(mean(sub$estimate), 3), ")")
         ),
         col = c("red", "orange", "darkgreen"),
         lty = 2, lwd = 2, cex = 0.7)
}

dev.off()

# ---- B: Density overlay ----
pdf("~/Desktop/3383/0_git/tiyani-aldex3/Results/effect_size_density_overlay.pdf", width=10, height=7)

comps <- unique(pvals$comparison)
colors <- rainbow(length(comps))

plot(NULL, xlim = range(pvals$estimate), ylim = c(0, 1),
     xlab = "Estimate (effect size)", ylab = "Density",
     main = "Effect Size Distributions Across All Pairwise Comparisons")

for (j in seq_along(comps)) {
  sub <- pvals[pvals$comparison == comps[j], ]
  lines(density(sub$estimate), col = colors[j], lwd = 1.5)
}

abline(v = 0, col = "black", lty = 2)
legend("topright", legend = comps, col = colors, lty = 1, lwd = 1.5, cex = 0.45)

dev.off()

# ---- C: Combined histogram ----
pdf("~/Desktop/3383/0_git/tiyani-aldex3/Results/effect_size_combined.pdf", width=8, height=6)

hist(pvals$estimate, breaks = 100,
     main = "Distribution of Effect Sizes Across All Pairwise Comparisons",
     xlab = "Estimate (effect size)",
     col = "steelblue", border = "white")
abline(v = 0, col = "red", lty = 2, lwd = 2)
abline(v = median(pvals$estimate), col = "orange", lty = 2, lwd = 2)
abline(v = mean(pvals$estimate), col = "darkgreen", lty = 2, lwd = 2)

legend("topright",
       legend = c(
         "Zero",
         paste0("Median (", round(median(pvals$estimate), 3), ")"),
         paste0("Mean (", round(mean(pvals$estimate), 3), ")")
       ),
       col = c("red", "orange", "darkgreen"),
       lty = 2, lwd = 2, cex = 0.7)

dev.off()

print("All three plots saved!")