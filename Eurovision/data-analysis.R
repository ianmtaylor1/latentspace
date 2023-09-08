# Code to do that data analysis on the Eurovision data

#library(amenhs)
library(coda)
library(ggplot2)
library(ggrepel)
library(foreach)
library(here)
library(knitr)
library(kableExtra)
library(dplyr)
library(tidyr)
library(reshape)

plotdir <- here("Eurovision", "plots")
resultdir <- here("Eurovision", "2015results")

# Read in all the csv files with response and covariates
Y <- read.csv(here("Eurovision", "2015data", "votes-ranks.csv"), row.names=1)
Xc.gender <- read.csv(here("Eurovision", "2015data", "covariate-gender.csv"), row.names=1, stringsAsFactors=FALSE)
Xd.lang <- read.csv(here("Eurovision", "2015data", "covariate-language.csv"), row.names=1)
Xd.contig <- read.csv(here("Eurovision", "2015data", "covariate-contig.csv"), row.names=1)
Xc.odds <- read.csv(here("Eurovision", "2015data", "covariate-odds.csv"), row.names=1)
Xc.pop <- read.csv(here("Eurovision", "2015data", "covariate-population.csv"), row.names=1)
Xc.gdp <- read.csv(here("Eurovision", "2015data", "covariate-gdp.csv"), row.names=1)

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

if (file.exists(file.path(resultdir, "no-rnd-effects.RDS"))) {
  res.no.re <- readRDS(file.path(resultdir, "no-rnd-effects.RDS"))
} else {
  res.no.re <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                           Xcol=Xc, Xdyad=Xd, 
                           family="rrl", rvar=FALSE, cvar=FALSE, dcor=FALSE,
                           plot=FALSE, gof=FALSE, print=TRUE,
                           nscan=nscan, burn=burn, odens=odens)
  saveRDS(res.no.re, file.path(resultdir, "no-rnd-effects.RDS"))
}

if (file.exists(file.path(resultdir, "standard.RDS"))) {
  res <- readRDS(file.path(resultdir, "standard.RDS"))
} else {
  res <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                     Xcol=Xc, Xdyad=Xd, 
                     family="rrl", rvar=FALSE, cvar=TRUE, dcor=FALSE,
                     plot=FALSE, gof=FALSE, print=TRUE,
                     nscan=nscan, burn=burn, odens=odens)
  saveRDS(res, file.path(resultdir, "standard.RDS"))
}

if (file.exists(file.path(resultdir, "hc-and-projection.RDS"))) {
  res.proj <- readRDS(file.path(resultdir, "hc-and-projection.RDS"))
} else {
  res.proj <- amenhs::ame(Y=as.matrix(Y[countries.in.final,countries.in.final]), 
                          Xcol=Xc, Xdyad=Xd,
                          family="rrl", rvar=FALSE, cvar=TRUE, halfcauchy=TRUE, project=TRUE, dcor=FALSE,
                          plot=FALSE, gof=FALSE, print=TRUE,
                          nscan=nscan, burn=burn, odens=odens)
  saveRDS(res.proj, file.path(resultdir, "hc-and-projection.RDS"))
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

ci.df <- foreach(i=1:length(covars), .combine="rbind") %do% {
  covar.name <- covars[i]
  display.name <- displaynames[i]
  
  covar.samples <- data.frame(
    Projected=c(rep("Restricted Network Model", nrow(res.proj$DELTA)), 
                rep("Non-Restricted Network Model", nrow(res.proj$BETA)),
                rep("No Random Effects", nrow(res.no.re$BETA))),
    Samples=c(res.proj$DELTA[,covar.name], 
              res.proj$BETA[,covar.name],
              res.no.re$BETA[,covar.name]),
    stringsAsFactors=FALSE)
  covar.ci <- data.frame(
    Covariate=rep(display.name, 3),
    Projected=c("Restricted Network Model","Non-Restricted Network Model", "No Random Effects"),
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
  
  #g <- ggplot(covar.samples, aes(x=Samples, fill=Projected)) + 
  #  geom_density(alpha=0.25) +
  #  ggtitle(display.name, subtitle="Posterior KDEs")
  #print(g)
  #g <- ggplot(covar.ci, aes(x=Projected, y=Mean, ymin=Low, ymax=High)) +
  #  geom_errorbar(width=0.2) +
  #  geom_point(size=3) +
  #  coord_cartesian(ylim=c(-1, 1.5)) +
  #  geom_abline(slope=0, intercept=0, color="red") +
  #  ggtitle(display.name, subtitle="Posterior means and 90% credible intervals")
  #print(g)
  covar.ci
}

ci.df$Projected <- factor(ci.df$Projected, levels=c("No Random Effects", "Non-Restricted Network Model", "Restricted Network Model"))
ci.df$Covariate <- factor(ci.df$Covariate, levels=c("Country Contiguity", "Log Betting Odds", "Log Population", "Log GDP per Capita"))

png(file.path(resultdir, "Eurovision-results-CI.png"), width=900, height=600)

ggplot(ci.df, aes(x=Covariate, y=Mean, ymin=Low, ymax=High, color=Projected, shape=Projected)) +
  geom_pointrange(position=position_dodge(width=0.5), size=1.2, fatten=3) +
  #geom_errorbar(size=1.5, width=0.3, position=position_dodge(width=0.5)) +
  #geom_point(size=5, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0) +
  #ggtitle("Posterior means and 90% credible intervals", 
  #        subtitle="Model Fixed Effects") +
  labs(y="Estimate", color="", shape="", x="") +
  theme_bw(base_family="serif") + 
  theme(text=element_text(size=25), legend.position=c(0.75, 0.85)) +
  coord_cartesian(ylim=c(-1.5, 1.5))

dev.off()


################################################################################
# Compare magnitude of posterior mean and CI width in table format
dp <- function(x, places) {
  format(round(x, places), nsmall=places)
}

ci.df <- ci.df |> mutate(Width=High - Low)

projected.ci.df <- ci.df |> filter(Projected=="Restricted Network Model")
other.ci.df <- ci.df |> filter(Projected!="Restricted Network Model")

projected.ci.df |> 
  inner_join(other.ci.df, by="Covariate") |>
  mutate(MeanRatio = Mean.x / Mean.y,
         WidthRatio = Width.x / Width.y,
         TableText = paste0(dp(MeanRatio, 3), ", ", dp(WidthRatio, 3))) |>
  select(Covariate, Model=Projected.y, MeanRatio, WidthRatio, TableText) |>
  pivot_wider(id_cols="Covariate", names_from="Model", values_from=c("TableText"), names_sort = TRUE) |>
  sort_df() |>
  kable(format="latex", booktabs=TRUE) |>
  add_header_above(c(" "=1,"Comparison Model"=2))

# Take 2 with separate columns for mean and width ratios
projected.ci.df |> 
  inner_join(other.ci.df, by="Covariate") |>
  mutate(MeanRatio = dp(Mean.x / Mean.y, 3),
         WidthRatio = dp(Width.x / Width.y, 3)) |>
  select(Covariate, Model=Projected.y, MeanRatio, WidthRatio) |>
  pivot_wider(id_cols="Covariate", names_from="Model", values_from=c("MeanRatio", "WidthRatio"), names_glue="{Model} {.value}", names_sort = TRUE) |>
  select(1,2,4,3,5) |> # Sort columns in the desired order
  kable(format="latex", booktabs=TRUE, 
        col.names=c("Covariate", "Mean Ratio", "Width Ratio", "Mean Ratio", "Width Ratio")) |>
  add_header_above(c(" "=1, "No Random Effects"=2, "Non-Restricted Network Model"=2)) |>
  add_header_above(c(" "=1,"Comparison Model"=4))


################################################################################
# Plot changes in random effect values (posterior means)

# pdf(file.path(resultdir, "Eurovision-results-plots.pdf"), width=8, height=8)
# 
# # Posterior means for column random effects with/without projections
# BPM <- data.frame(
#   Projected=res.proj$BPM.ORTH,
#   NotProjected=res.proj$BPM,
#   Country=names(res.proj$BPM),
#   stringsAsFactors=FALSE
# )
# BPM$AbsChg <- abs(BPM$Projected - BPM$NotProjected)
# BPM$ChgRnk <- rank(BPM$AbsChg)
# BPM$TopMovers <- ifelse(BPM$ChgRnk >= 23, BPM$Country, NA)
# 
# # Posterior mean of column effect scatterplots
# ggplot(BPM, aes(x=NotProjected, y=Projected)) + 
#   geom_point() +
#   ggtitle("Column Random Effects", subtitle="Posterior means before and after projections") +
#   geom_text_repel(aes(label=TopMovers), size=3.5) +
#   geom_abline(slope=1, intercept=0, color="red")
# # Absolute change by country
# ggplot(BPM, aes(x=Country, y=AbsChg)) + 
#   geom_bar(stat="identity")
# 
# dev.off()