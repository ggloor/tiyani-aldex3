library(ALDEx2)
library(ALDEx3)

# Load all datasets
load('~/Desktop/3383/0_git/tiyani-aldex3/data/cent.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/eld.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/kin.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mage.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/mid.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/pup.age.Rda')
load('~/Desktop/3383/0_git/tiyani-aldex3/data/you.age.Rda')




Y <- cbind(kin,pup)
# make a conditions vector
# for a continuous variable (age) we replace this with the ages of the samples
age <- c(as.numeric(kin.age[1,]), as.numeric(pup.age[1,]))
conds <- c(rep("K", 103), rep('P',161))
data <- data.frame(age=age, condition=conds)

print("Test_block_2")

# does the calculation
#  ## Had to add "ALDEx3::aldex()" to make the code work: for some reason it was not able to work 
#  ## and started giving me the "Error in round(conds) : non-numeric argument to mathematical function"
#  ## Till I made that change
res.age <- ALDEx3::aldex(Y, ~age, data, nsample=128, scale=clr.sm, gamma=0.3)
res.conds <- ALDEx3::aldex(Y, ~condition, data, nsample=128, scale=clr.sm, gamma=1e-3)

# summarize
# The command was giving a error "could not find the function" so, It has been updated to newer command 
# summary.aldex() --> summary()
sum.0 <- summary(res.conds)
sum.0a <- summary(res.age)
#plot(sum.0$estimate, sum.0a$estimate)

#Checkpoint to validate the pervious code is working
print("Test_Block_3")

sig <- which(sum.0$p.val.adj < 0.05)
siga <- which(sum.0a$p.val.adj < 0.05)

plot(sum.0a$std.error, sum.0a$estimate)

points(sum.0a$std.error[sig], sum.0a$estimate[sig], col='red', pch=19, cex=0.5)
points(sum.0a$std.error[siga], sum.0a$estimate[siga], col='orange', pch=19, cex=0.5)
abline(h=0)
