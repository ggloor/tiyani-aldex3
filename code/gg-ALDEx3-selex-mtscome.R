# this adds a constant to the TSS scale model
const.sm <- function(X, logComp, const=0, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  LambdaScale <- matrix(rnorm(P*nsample,const,gamma), P, nsample)
  logScale <- t(X)%*% LambdaScale
  return(logScale)
}

selex <- read.table('~/Documents/0_git/datasets/selex.txt', header=T, 
	row.names=1, sep='\t')

# estimate G from log2 proportions	
apply(selex+0.5, 2, function(x) mean(log2(x/sum(x))))
# suggests ~2e7 difference in scale between groups

# estimate mean in log2
# equivalent to 1/D
apply(selex+0.5, 2, function(x) log2(mean(x/sum(x))))
# suggests each group has equivalent scale

Y <- selex
conds <- c(rep("N", 7), rep("S",7))
X <- data.frame(conds))

sel.clr <- aldex(Y, ~conds, X, nsample=128, scale=clr.sm, gamma=0)
sel.clr2 <- aldex(Y, ~conds, X, nsample=128, scale=clr.sm, gamma=0)
aldex.plot(sel.clr, contrast='conds', plot='effect')
clr.eff <- aldex.effect(sel.clr, contrast='conds')
clr.eff2 <- aldex.effect(sel.clr2, contrast='conds')

sel.tss <- aldex(Y, ~conds, X, nsample=128, scale=const.sm, gamma=0)
tss.eff <- aldex.effect(sel.tss, contrast='conds')

# hist of clr.eff vs tss.eff, suggests modal diff of -8 for tss
# so add const=8 and re-run
sel.tss8 <- aldex(Y, ~conds, X, nsample=128, scale=const.sm, const=8, gamma=0)
sel.tss0 <- aldex(Y, ~conds, X, nsample=128, scale=const.sm, const=0, gamma=0)
tss.eff8 <- aldex.effect(sel.tss8, contrast='conds')
tss.eff0 <- aldex.effect(sel.tss0, contrast='conds')


##### ALDEx3 on vaginal dataset
### WARNING: paths in setup.R need changing along with hard-coded paths below
### REQUIRES ALDEx3 which can be obtained from CRAN
library(ALDEx3)

devtools::load_all('~/Documents/0_git/CoDaSeq/CoDaSeq')
path.to.github <- "~/Documents/0_git/projects/dossantos2024study/"

# set path in setup.R to gg
source(paste(path.to.github, "code/setup.R", sep = ""))

# load in vNumber -> KEGG pathway lookup table from VIRGO
path.table <- read.table(paste(locn,'1_VIRGO/8.C.kegg.pathway.copy.txt', sep=""), 
                         sep="\t", header=T, row.names=1, fill=TRUE)

# load in vector containing cst info from london/europe species heatmap
load(paste(path.to.github, 'Rdata/hm.metadata.Rda', sep = ""))

# load in london/europe heatmap column colour bar list
load(paste(path.to.github, "Rdata/hm.column.cols.Rda",sep = ""))

# remove two samples from the filtered london/europe feature table aggregated by
# K0 number (classed as BV but almost no BV organisms):
#   -  v.001A: close to 100% L. gasseri with practically no BV organisms
#   -  v.019A: around 80 % iners, ~ 20% crispatus and a tiny bit of Gardnerella
ko.both<-ko.both[,-c(9,22)]

# remove non-bacterial KOs from the K number-aggregated, filtered feature table
# (these were discovered during curation of 'Unknown' pathways)
#    - K03364: Eukaryotic cell division cycle 20-like protein 1
#    - K13963: Serpin B (eukaryotic serine protease inhibitor)
#    - K01173: Mitochondrial endonuclease G
#    - K12373: Human lysosomal hexosaminidase
#    - K14327: Eukaryotic regulator of nonsense transcripts 2
#    - K00863: Human triose/dihydroxyacetone kinase
#    - K00599: Eukaryotic tRNA N(3)-methylcytidine methyltransferase
#    - K13993: Human HSP20
#    - K00811: Chloroplastic aspartate aminotransferase
#    - K03260: Eukaryotic translation initiation factor 4G
#    - K00985: Enterovirus RNA-directed RNA polymerase

ko.both <- ko.both[-which(grepl(paste("K03364","K13963","K01173","K12373",
                                      "K14327","K00863","K00599","K13993",
                                      "K00811","K03260","K00985", sep = "|"),
                                rownames(ko.both))),]
# clr.sm modified to accept a constant for the offset of the estimate
clr.const.sm <- function(X, logComp, const=0, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  logScale <- -colMeans(logComp, dims=1)

  tmp <- P*nsample
  LambdaScale <- matrix(rnorm(tmp,const,gamma), P, nsample)
  logScale <- logScale + t(X)%*% LambdaScale
  return(logScale)
}


ko.conds <- c(rep('H',8), rep('B',12), rep('B',14), rep('H', 8)) 

X <- data.frame(ko.conds)

# generate ALDEx3 outputs for plotting
vag.clr <- aldex(ko.both, ~ko.conds, X, nsample=128, scale=clr.sm, gamma=0)
vag.clr5 <- aldex(ko.both, ~ko.conds, X, nsample=128, scale=clr.const.sm, const=-3, gamma=0.5)
vag.tss <- aldex(ko.both, ~ko.conds, X, nsample=128, scale=tss.sm, gamma=0)
vag.tss5 <- aldex(ko.both, ~ko.conds, X, nsample=128, scale=tss.sm, gamma=0.5)

# shows that tss is the same as clr with offset of -3
par(mfrow=c(2,2))
aldex.plot(vag.clr, contrast='ko.conds', main='clr g=0', plot='eff', cex=0.6)
#abline(h=3, lty=2, col='grey')
aldex.plot(vag.tss, contrast='ko.conds', main='tss g=0', plot='eff', cex=0.6)
aldex.plot(vag.clr5, contrast='ko.conds', main='clr c=-3, g=0.5', plot='eff', cex=0.6)
aldex.plot(vag.tss5, contrast='ko.conds', main='tss g=0.5', plot='eff', cex=0.6)

sum.clr5 <- summary(vag.clr5)
sum.tss5 <- summary(vag.tss5)



