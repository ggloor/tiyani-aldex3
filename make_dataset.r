# NOTE
# this dataset included only samples from age 3 onward, and excludes duplicates and
# extraction controls, etc
#
# set a few parameters before calling this.
# addys TRUE or FALSE # determines if ys are included in the ouput or notinf.clr
 min.prop = 0.001
 min.occurrence = 0.0

# do not include young soldier group
addys <- TRUE

# load the required R packages
require(compositions) # exploratory data analysis of compositional data
require(zCompositions) # used for 0 substitution
require(ALDEx2) # used for per-OTU comparisons
source('code/codaSeq.filter.r' )
require(xtable) # used to generate tables from datasets
library(igraph) # used to generate graphs from phi data
library(car) # used to generate graphs from phi data
#library(propr) # used for correlation calculation
library(vegan)
# you will need to download this directly from github
# I put the files in a directory called git at the base of my user space
# https://github.com/ggloor/propr
# eventually, these will be in codaSeq package
# source("~/git/CoDaSeq/chunk/codaSeq_functions.R")
# source("~/git/propr/R/propr-functions.R")
# source("chunk/separability_measures.R")

colours <- c("indianred1", "steelblue3",  "skyblue1", "mediumorchid","royalblue4", "olivedrab3",
   "pink", "#FFED6F", "mediumorchid4", "ivory2", "tan1", "aquamarine3", "#C0C0C0",
    "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF",
   "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC",
   "#669999", "#CCCC66", "#CC6600", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00",
   "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED"
)

rbcol <- rainbow(n=10, alpha=0.5)

# load in the metadata
m <- read.table("~/Documents/0_git/tianyi/data/clean_meta.txt", header=T, row.names=1, sep="\t",
    check.names=F, comment.char="",quote="")

# load in the table of OTU counts
d <- read.table("~/Documents/0_git/tianyi/data/clean_data.txt",header=T, row.names=1, sep="\t",
    check.names=F, comment.char="",quote="", stringsAsFactors=F)

# load in taxonomy information
t.df <- read.table("~/Documents/0_git/tianyi/data/clean_tax.txt",header=T, row.names=1, sep="\t",
    check.names=F, comment.char="",quote="", stringsAsFactors=F)

d.m <- d[,colnames(d) %in% rownames(m)]

# set dataset to include or exclude ys cohort
# kindergarten, Pupils, mid_school, youth, young soldier, mid_age, elder, Centenarians
# 3:6, 8:12, 13:14, 19:22, 19:29, 30:50, 60:80, >90


kin <- d.m[, rownames(m)[as.numeric(m$Age) > 2 & as.numeric(m$Age) < 7 & m$Group == 'kindergarten']]
pup <- d.m[, rownames(m)[as.numeric(m$Age) > 7 & as.numeric(m$Age) < 13 & m$Group == 'Pupils']]
pup.age <-  data.frame(rbind(m[colnames(pup), "Age"]))
colnames(pup.age) <- colnames(pup)

mid <- d.m[, rownames(m)[as.numeric(m$Age) > 12 & as.numeric(m$Age) < 15 & m$Group == 'mid_school']]
you <- d.m[, rownames(m)[as.numeric(m$Age) > 18 & as.numeric(m$Age) < 25 & m$Group == 'youth']]
if(addys == TRUE) ys <- d.m[, rownames(m)[as.numeric(m$Age) > 18 & as.numeric(m$Age) < 25 & m$Group == 'young soldier']]
mage <- d.m[, colnames(d.m) %in% rownames (m)[as.numeric(m$Age) > 29 & as.numeric(m$Age) < 51 & m$Group == 'mid_age']]
eld <- d.m[, colnames(d.m) %in% rownames(m)[as.numeric(m$Age) > 59 & as.numeric(m$Age) < 80 & m$Group == 'elder']]
cent <- d.m[, colnames(d.m) %in% rownames(m)[as.numeric(m$Age) > 94 & m$Group == 'Centenarians']]

if(addys == TRUE) clean.data <- data.frame(kin,pup,mid,you,ys,mage,eld,cent, stringsAsFactors=F)
if(addys == FALSE) clean.data <- data.frame(kin,pup,mid,you,mage,eld,cent, stringsAsFactors=F)

if(addys == TRUE) legend.names = c("3-6", "8-12", "13-14", "19-24", "19-24S", "30-50", "60-79", ">94")
if(addys == FALSE) legend.names = c("3-6", "8-12", "13-14", "19-24", "30-50", "60-79", ">94")


# remove 0 count OTUs
d.notinf <-  codaSeq.filter(clean.data, min.prop=0, min.occurrence=0, samples.by.row=FALSE)
colnames(d.notinf) <- gsub("^X", "", colnames(d.notinf))

# keep only OTUs that are >0.001 abundant across all samples 
# output is samples by column. 1391 OTUs, 1095 samples
notinf.filt <- codaSeq.filter(d.notinf,min.prop=min.prop, max.prop=1, min.occurrence=min.occurrence	, samples.by.row=FALSE)

# make the data table
kin <- kin[rownames(notinf.filt),]
# pull the ages
kin.age <-  data.frame(rbind(m[colnames(kin), "Age"]))
# rename both
colnames(kin) <- paste('K_', colnames(kin), sep="")
colnames(kin.age) <- colnames(kin)

pup <- pup[rownames(notinf.filt),]
pup.age <-  data.frame(rbind(m[colnames(pup), "Age"]))
colnames(pup) <- paste('P_', colnames(pup), sep="")
colnames(pup.age) <- colnames(pup)

mid <- mid[rownames(notinf.filt),]
mid.age <-  data.frame(rbind(m[colnames(mid), "Age"]))
colnames(mid) <- paste('M_', colnames(mid), sep="")
colnames(mid.age) <- colnames(mid)

you <- you[rownames(notinf.filt),]
you.age <-  data.frame(rbind(m[colnames(you), "Age"]))
colnames(you) <- paste('Y_', colnames(you), sep="")
colnames(you.age) <- colnames(you)

if(addys == TRUE) {
  ys <- ys[rownames(notinf.filt),]
  ys.age <-  data.frame(rbind(m[colnames(ys), "Age"]))
  colnames(ys) <- paste('S_', colnames(ys), sep="")
  colnames(ys.age) <- colnames(ys)

}
mage <- mage[rownames(notinf.filt),]
mage.age <-  data.frame(rbind(m[colnames(mage), "Age"]))
colnames(mage) <- paste('A_', colnames(mage), sep="")
colnames(mage.age) <- colnames(mage)

eld <- eld[rownames(notinf.filt),]
eld.age <-  data.frame(rbind(m[colnames(eld), "Age"]))
colnames(eld) <- paste('E_', colnames(eld), sep="")
colnames(eld.age) <- colnames(eld)

cent <- cent[rownames(notinf.filt),]
cent.age <-  data.frame(rbind(m[colnames(cent), "Age"]))
colnames(cent) <- paste('C_', colnames(cent), sep="")
colnames(cent.age) <- colnames(cent)

save(kin, file="data/kin.Rda")
save(pup, file="data/pup.Rda")
save(mid, file="data/mid.Rda")
save(you, file="data/you.Rda")
save(mage, file="data/mage.Rda")
save(eld, file="data/eld.Rda")
save(cent, file="data/cent.Rda")

save(kin.age, file="data/kin.age.Rda")
save(pup.age, file="data/pup.age.Rda")
save(mid.age, file="data/mid.age.Rda")
save(you.age, file="data/you.age.Rda")
save(mage.age, file="data/mage.age.Rda")
save(eld.age, file="data/eld.age.Rda")
save(cent.age, file="data/cent.age.Rda")


