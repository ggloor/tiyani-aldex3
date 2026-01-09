# ALDEx2

load('data/kin.Rda')
load("data/pup.Rda")

library(ALDEx2)
conds <- c(rep("k",103), rep("P",161))
KP.0 <- aldex(cbind(kin,pup), conditions=conds, gamma=0.2)
aldex.plot(KP.0) # effect plot
aldex.plot(KP.0, type='volcano')

# ALDEx3
# loading a non-compiled version from local github repository
devtools::load_all('~/Documents/0_git/projects/ALDEx3')

library(ALDEx3)
#      Y <- matrix(1:110, 10, 11)
#      condition <- c(rep(0, 5), rep(1, 6))
#      data <- data.frame(condition=condition)
#      ## demonstrate formula interface and passing optional argument (gamma) to
#      ## the scale model (clr)
#      res <- aldex(Y, ~condition, data, nsample=2000, scale=clr.sm, gamma=0.5)

# make the data matrix via cbind
Y <- cbind(kin,pup)
# make a conditions vector
# for a continuous variable (age) we replace this with the ages of the samples
conds <- c(rep("K", 103), rep('P',161))
data <- data.frame(condition=conds)


# does the calculation
res <- aldex(Y, ~condition, data, nsample=128, scale=clr.sm, gamma=1e-3)

# summarize
sum.0 <- summary.aldex(res)

# there is not a native plotting function yet!!
plot(sum.0$std.error, sum.0$estimate)
sig <- which(sum.0$p.val.adj < 0.05)
points(sum.0$std.error[sig], sum.0$estimate[sig], col='red', pch=19, cex=0.5)
abline(h=0)

# now do for gamma of 0.3
res <- aldex(Y, ~condition, data, nsample=128, scale=clr.sm, gamma=0.3)

sum.3 <- summary.aldex(res)
plot(sum.3$std.error, sum.3$estimate, col="grey", pch=19, cex=0.5)
sig3 <- which(sum.3$p.val.adj < 0.05)
# these are signficant without gamma
points(sum.3$std.error[sig], sum.3$estimate[sig], col='orange', pch=19, cex=0.5)
# these are signficant with gamma
points(sum.3$std.error[sig3], sum.3$estimate[sig3], col='red', pch=19, cex=0.5)
abline(h=0)





