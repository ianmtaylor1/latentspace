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

if (file.exists("2015results/no-rnd-effects.RDS")) {
  res.no.re <- readRDS("2015results/no-rnd-effects.RDS")
} else {
  res.no.re <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                           Xcol=Xc, Xdyad=Xd, 
                           family="rrl", rvar=FALSE, cvar=FALSE, dcor=FALSE,
                           plot=FALSE, gof=FALSE, print=TRUE,
                           nscan=500000, burn=10000, odens=500)
  saveRDS(res.no.re, "2015results/no-rnd-effects.RDS")
}

if (file.exists("2015results/standard.RDS")) {
  res <- readRDS("2015results/standard.RDS")
} else {
  res <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                     Xcol=Xc, Xdyad=Xd, 
                     family="rrl", rvar=FALSE, cvar=TRUE, dcor=FALSE,
                     plot=FALSE, gof=FALSE, print=TRUE,
                     nscan=500000, burn=10000, odens=500)
  saveRDS(res, "2015results/standard.RDS")
}

if (file.exists("2015results/hc-and-projection.RDS")) {
  res.proj <- readRDS("2015results/hc-and-projection.RDS")
} else {
  res.proj <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                          Xcol=Xc, Xdyad=Xd,
                          family="rrl", rvar=FALSE, cvar=TRUE, halfcauchy=TRUE, project=TRUE, dcor=FALSE,
                          plot=FALSE, gof=FALSE, print=TRUE,
                          nscan=500000, burn=10000, odens=500)
  saveRDS(res.proj, "2015results/hc-and-projection.RDS")
}


# Summarize results (BETA)
cat("No Random Effects\n\n")
colMeans(res.no.re$BETA)
effectiveSize(res.no.re$BETA)
apply(res.no.re$BETA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

cat("Standard Random Effects\n\n")
colMeans(res$BETA)
effectiveSize(res$BETA)
apply(res$BETA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

cat("With Half-Cauchy\n\n")
colMeans(res.proj$BETA)
effectiveSize(res.proj$BETA)
apply(res.proj$BETA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

cat("With HC and Projected\n\n")
colMeans(res.proj$DELTA)
effectiveSize(res.proj$DELTA)
apply(res.proj$DELTA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
