# Code to do that data analysis on the Eurovision data

library(amenhs)
library(coda)

# Read in all the csv files with response and covariates
Y <- read.csv("2015data/votes-ranks.csv", row.names=1)
Xc.gender <- read.csv("2015data/covariate-gender.csv", row.names=1, stringsAsFactors=FALSE)
Xd.lang <- read.csv("2015data/covariate-language.csv", row.names=1)
Xd.contig <- read.csv("2015data/covariate-contig.csv", row.names=1)
Xc.odds <- read.csv("2015data/covariate-odds.csv", row.names=1)

# Which countries are we restricting the analysis to?
countries.in.final <- row.names(Xc.gender)[Xc.gender[,"Female"] != "#N/A"]

Xd <- array(
  c(#as.matrix(Xd.lang[countries.in.final, countries.in.final]),
    as.matrix(Xd.contig[countries.in.final, countries.in.final])),
  dim=c(length(countries.in.final), length(countries.in.final), 1)
)
dimnames(Xd)[[3]] <- list("Contig")
#Xc <- data.matrix(cbind(
#  Xc.gender[countries.in.final,],
#  Xc.odds[countries.in.final,"Median16"]
#))
#colnames(Xc) <- c(colnames(Xc.gender), "MedianOdds")
Xc <- data.matrix(data.frame(MedianOdds=Xc.odds[countries.in.final,"Median16"]))

res.no.re <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                   Xcol=Xc, Xdyad=Xd, 
                   family="rrl", rvar=FALSE, cvar=FALSE, dcor=FALSE,
                   plot=FALSE, gof=FALSE, print=TRUE,
                   nscan=500000, burn=10000, odens=500)

res <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                   Xcol=Xc, Xdyad=Xd, 
                   family="rrl", rvar=FALSE, cvar=TRUE, dcor=FALSE,
                   plot=FALSE, gof=FALSE, print=TRUE,
                   nscan=500000, burn=10000, odens=500)

#res.hc <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
#                   Xcol=Xc, Xdyad=Xd, 
#                   family="rrl", rvar=FALSE, cvar=TRUE, halfcauchy=TRUE, dcor=FALSE,
#                   plot=FALSE, gof=FALSE, print=FALSE,
#                   nscan=500000, burn=10000, odens=500)

res.proj <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                   Xcol=Xc, Xdyad=Xd,
                   family="rrl", rvar=FALSE, cvar=TRUE, halfcauchy=TRUE, project=TRUE, dcor=FALSE,
                   plot=FALSE, gof=FALSE, print=TRUE,
                   nscan=500000, burn=10000, odens=500)

# Summarize results (BETA)
cat("No Random Effects\n")
colMeans(res.no.re$BETA)
effectiveSize(res.no.re$BETA)
apply(res.no.re$BETA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

cat("Standard Random Effects\n")
colMeans(res$BETA)
effectiveSize(res$BETA)
apply(res$BETA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

cat("With Half-Cauchy")
colMeans(res.proj$BETA)
effectiveSize(res.proj$BETA)
apply(res.proj$BETA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

cat("With HC and Projected")
colMeans(res.proj$DELTA)
effectiveSize(res.proj$DELTA)
apply(res.proj$DELTA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
