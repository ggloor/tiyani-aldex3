# ALDEx2

load('~/Desktop/3383/0_git/tiyani-aldex3/data/kin.Rda')
load("~/Desktop/3383/0_git/tiyani-aldex3/data/pup.Rda")

library(ALDEx2)
conds <- c(rep("k",103), rep("P",161))
#KP.0 <- aldex(cbind(kin,pup), conditions=conds, gamma=0.2)
#aldex.plot(KP.0) # effect plot
#aldex.plot(KP.0, type='volcano')

print("ALDEX2-done") #Test output

# ALDEx3
# loading a non-compiled version from local github repository
#devtools::load_all('~/Documents/0_git/projects/ALDEx3')

library(ALDEx3)
      Y <- matrix(1:110, 10, 11)
      condition <- c(rep(0, 5), rep(1, 6))
      data <- data.frame(condition=condition)
#      ## demonstrate formula interface and passing optional argument (gamma) to
#      ## the scale model (clr)
      
      res <- ALDEx3::aldex(Y, ~condition, data, nsample=2000, scale=clr.sm, gamma=0.5)
#  ## Had to add "ALDEx3::aldex()" to make the code work: for some reason it was not able to work 
#  ## and started giving me the "Error in round(conds) : non-numeric argument to mathematical function"
#  ## Till I made that change

print("Test_Block_1")


# make the data matrix via cbind
Y <- cbind(kin,pup)
# make a conditions vector
# for a continuous variable (age) we replace this with the ages of the samples
conds <- c(rep("K", 103), rep('P',161))
data <- data.frame(condition=conds)

print("Test_block_2")

# does the calculation
#  ## Had to add "ALDEx3::aldex()" to make the code work: for some reason it was not able to work 
#  ## and started giving me the "Error in round(conds) : non-numeric argument to mathematical function"
#  ## Till I made that change
res <- ALDEx3::aldex(Y, ~condition, data, nsample=128, scale=clr.sm, gamma=1e-3)

# summarize
# The command was giving a error "could not find the function" so, It has been updated to newer command 
# summary.aldex() --> summary()
sum.0 <- summary(res)

#Checkpoint to validate the pervious code is working
print("Test_Block_3")

# there is not a native plotting function yet!!
plot(sum.0$std.error, sum.0$estimate)
sig <- which(sum.0$p.val.adj < 0.05)
points(sum.0$std.error[sig], sum.0$estimate[sig], col='red', pch=19, cex=0.5)
abline(h=0)

print("Test_block_4")

## now do for gamma of 0.3
#res <- ALDEx3::aldex(Y, ~condition, data, nsample=128, scale=clr.sm, gamma=0.3)
#  ## Had to add "ALDEx3::aldex()" to make the code work: for some reason it was not able to work 
#  ## and started giving me the "Error in round(conds) : non-numeric argument to mathematical function"
#  ## Till I made that change


## The command was giving a error "could not find the function" so, It has been updated to newer command 
## summary.aldex() --> summary()
#sum.3 <- summary(res)
#plot(sum.3$std.error, sum.3$estimate, col="grey", pch=19, cex=0.5)
#sig3 <- which(sum.3$p.val.adj < 0.05)

## these are signficant without gamma
#points(sum.3$std.error[sig], sum.3$estimate[sig], col='orange', pch=19, cex=0.5)

## these are signficant with gamma
#points(sum.3$std.error[sig3], sum.3$estimate[sig3], col='red', pch=19, cex=0.5)
#abline(h=0)



