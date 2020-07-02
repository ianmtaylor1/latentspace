# Code to do that data analysis on the Eurovision data

library(amenhs)

# Read in all the csv files with response and covariates
Y <- read.csv("2015data/votes-ranks.csv", row.names=1)
Xc.gender <- read.csv("2015data/covariate-gender.csv", row.names=1, stringsAsFactors=FALSE)
Xd.lang <- read.csv("2015data/covariate-language.csv", row.names=1)

# Which countries are we restricting the analysis to?
countries.in.final <- row.names(Xc.gender)[Xc.gender[,"Female"] != "#N/A"]

res.no.re <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                   Xcol=data.matrix(Xc.gender[countries.in.final,]), 
                   Xdyad=as.array(as.matrix(Xd.lang[countries.in.final,countries.in.final])), 
                   family="frn", rvar=FALSE, cvar=FALSE)

res <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                   Xcol=data.matrix(Xc.gender[countries.in.final,]), 
                   Xdyad=as.array(as.matrix(Xd.lang[countries.in.final,countries.in.final])), 
                   family="frn")

res.hc <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                   Xcol=data.matrix(Xc.gender[countries.in.final,]), 
                   Xdyad=as.array(as.matrix(Xd.lang[countries.in.final,countries.in.final])), 
                   family="frn", halfcauchy=TRUE)

res.proj <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                   Xcol=data.matrix(Xc.gender[countries.in.final,]), 
                   Xdyad=as.array(as.matrix(Xd.lang[countries.in.final,countries.in.final])), 
                   family="frn", halfcauchy=TRUE, project=TRUE)

# Summarize results (BETA)
colMeans(res$BETA)
apply(res$BETA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
colMeans(res.hc$BETA)
apply(res.hc$BETA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
colMeans(res.proj$BETA)
apply(res.proj$BETA, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))