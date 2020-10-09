# Code to do that data analysis on the Eurovision data

library(amenhs)
library(coda)
library(ggplot2)
library(ggrepel)

# Read in all the csv files with response and covariates
Y <- read.csv("2015data/votes-ranks.csv", row.names=1)
Xc.gender <- read.csv("2015data/covariate-gender.csv", row.names=1, stringsAsFactors=FALSE)
Xd.lang <- read.csv("2015data/covariate-language.csv", row.names=1)
Xd.contig <- read.csv("2015data/covariate-contig.csv", row.names=1)
Xc.odds <- read.csv("2015data/covariate-odds.csv", row.names=1)
Xc.pop <- read.csv("2015data/covariate-population.csv", row.names=1)
Xc.gdp <- read.csv("2015data/covariate-gdp.csv", row.names=1)

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
Xc <- data.matrix(data.frame(LogMedianOdds=log(Xc.odds[countries.in.final,"Median16"]),
                             LogPopulation=log(Xc.pop[countries.in.final,"Pop2015"]),
                             LogGDP=log(Xc.gdp[countries.in.final,"GDPpc2015"])))

# Run specifics
nscan <- 5000000
burn <- 25000
odens <- 1000

if (file.exists("2015results/no-rnd-effects.RDS")) {
  res.no.re <- readRDS("2015results/no-rnd-effects.RDS")
} else {
  res.no.re <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                           Xcol=Xc, Xdyad=Xd, 
                           family="rrl", rvar=FALSE, cvar=FALSE, dcor=FALSE,
                           plot=FALSE, gof=FALSE, print=TRUE,
                           nscan=nscan, burn=burn, odens=odens)
  saveRDS(res.no.re, "2015results/no-rnd-effects.RDS")
}

if (file.exists("2015results/standard.RDS")) {
  res <- readRDS("2015results/standard.RDS")
} else {
  res <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                     Xcol=Xc, Xdyad=Xd, 
                     family="rrl", rvar=FALSE, cvar=TRUE, dcor=FALSE,
                     plot=FALSE, gof=FALSE, print=TRUE,
                     nscan=nscan, burn=burn, odens=odens)
  saveRDS(res, "2015results/standard.RDS")
}

if (file.exists("2015results/hc-and-projection.RDS")) {
  res.proj <- readRDS("2015results/hc-and-projection.RDS")
} else {
  res.proj <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                          Xcol=Xc, Xdyad=Xd,
                          family="rrl", rvar=FALSE, cvar=TRUE, halfcauchy=TRUE, project=TRUE, dcor=FALSE,
                          plot=FALSE, gof=FALSE, print=TRUE,
                          nscan=nscan, burn=burn, odens=odens)
  saveRDS(res.proj, "2015results/hc-and-projection.RDS")
}

################################################################################

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


################################################################################

# Make plots of CIs for each covariate
covars <- c("LogMedianOdds.col", "LogPopulation.col", "LogGDP.col", ".dyad")
displaynames <- c("Log Betting Odds", "Log Population", "Log GDP per Capita", "Country Contiguity")

pdf("2015results/Eurovision-results-plots.pdf", width=8, height=8)

for (i in 1:length(covars)) {
  covar.name <- covars[i]
  display.name <- displaynames[i]
  
  covar.samples <- data.frame(
    Projected=c(rep("Projected", nrow(res.proj$DELTA)), 
                rep("Not Projected", nrow(res.proj$BETA)),
                rep("No Random Effects", nrow(res.no.re$BETA))),
    Samples=c(res.proj$DELTA[,covar.name], 
              res.proj$BETA[,covar.name],
              res.no.re$BETA[,covar.name]),
    stringsAsFactors=FALSE)
  covar.ci <- data.frame(
    Projected=c("Projected","Not Projected", "No Random Effects"),
    Mean=c(mean(res.proj$DELTA[,covar.name]), 
           mean(res.proj$BETA[,covar.name]),
           mean(res.no.re$BETA[,covar.name])),
    Low=c(quantile(res.proj$DELTA[,covar.name], probs=0.05), 
          quantile(res.proj$BETA[,covar.name], probs=0.05),
          quantile(res.no.re$BETA[,covar.name], probs=0.05)),
    High=c(quantile(res.proj$DELTA[,covar.name], probs=0.95), 
           quantile(res.proj$BETA[,covar.name], probs=0.95),
           quantile(res.no.re$BETA[,covar.name], probs=0.95)),
    stringsAsFactors=FALSE
  )
  
  g <- ggplot(covar.samples, aes(x=Samples, fill=Projected)) + 
    geom_density(alpha=0.25) +
    ggtitle(display.name, subtitle="Posterior KDEs")
  print(g)
  g <- ggplot(covar.ci, aes(x=Projected, y=Mean, ymin=Low, ymax=High)) +
    geom_errorbar(width=0.2) +
    geom_point(size=1.5) +
    ggtitle(display.name, subtitle="Posterior means and 90% credible intervals")
  print(g)
}


# Posterior means for column random effects with/without projections
BPM <- data.frame(
  Projected=res.proj$BPM.ORTH,
  NotProjected=res.proj$BPM,
  Country=names(res.proj$BPM),
  stringsAsFactors=FALSE
)
BPM$AbsChg <- abs(BPM$Projected - BPM$NotProjected)
BPM$ChgRnk <- rank(BPM$AbsChg)
BPM$TopMovers <- ifelse(BPM$ChgRnk >= 23, BPM$Country, NA)

# Posterior mean of column effect scatterplots
ggplot(BPM, aes(x=NotProjected, y=Projected)) + 
  geom_point() +
  ggtitle("Column Random Effects", subtitle="Posterior means before and after projections") +
  geom_text_repel(aes(label=TopMovers), size=3.5) +
  geom_abline(slope=1, intercept=0, color="red")
# Absolute change by country
ggplot(BPM, aes(x=Country, y=AbsChg)) + 
  geom_bar(stat="identity")

dev.off()