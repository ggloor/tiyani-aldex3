load('data/kin.Rda')
load("data/pup.Rda")

library(ALDEx2)
conds <- c(rep("k",103), rep("P",161))
KP.0 <- aldex(cbind(kin,pup), conditions=conds, gamma=0.2)
aldex.plot(KP.0)

