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




