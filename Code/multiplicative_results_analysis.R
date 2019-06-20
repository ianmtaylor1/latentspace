# Analyze the results of the first Simulation using the custom version of the code
library(RColorBrewer)
library(ggplot2)

setwd("C:/Users/Ian Taylor/Documents/Git/latentspace")

infile <- "Results/Multiplicative Results n=20.csv"
outfile <- "Results/Multiplicative Results n=20.pdf"

mypalette <- c(
  brewer.pal(10,"Paired")[c(1,3,5,9)]
)

palette(mypalette)


results <- read.csv(infile, header=TRUE, stringsAsFactors=FALSE)

results[,"Covered"] <- 1 * ((results[,"TrueValue"] >= results[,"CI_low"]) & (results[,"TrueValue"] <= results[,"CI_high"]))
results[,"CI_width"] <- results[,"CI_high"] - results[,"CI_low"]

aggregate <- unique(results[,c("Run","Xseed","Z_additive","Z_multiplicative","AddRE","MulREdim","Variable","TrueValue")])

for (r in 1:(dim(aggregate)[1])) {
  subresults <- results[(results[,"Run"] == aggregate[r,"Run"]) # Xseed is based on Run
                        & (results[,"Z_additive"] == aggregate[r,"Z_additive"])
                        & (results[,"Z_multiplicative"] == aggregate[r,"Z_multiplicative"])
                        & (results[,"AddRE"] == aggregate[r,"AddRE"])
                        & (results[,"MulREdim"] == aggregate[r,"MulREdim"])
                        & (results[,"Variable"] == aggregate[r,"Variable"]) # TrueValue is based on Variable
                        ,]
  aggregate[r,"n_runs"] <- dim(subresults)[1]
  aggregate[r,"n_covered"] <- sum(subresults[,"Covered"])
  aggregate[r,"estimate_mean"] <- mean(subresults[,"Estimate"])
  aggregate[r,"estimate_var"] <- var(subresults[,"Estimate"])
  aggregate[r,"mean_CI_width"] <- mean(subresults[,"CI_width"])
  if((aggregate[r,"AddRE"] == FALSE) && (aggregate[r,"MulREdim"] == 0)) {
    aggregate[r,"RndEffects"] <- "None"
  } else if ((aggregate[r,"AddRE"] == TRUE) && (aggregate[r,"MulREdim"] == 0)) {
    aggregate[r,"RndEffects"] <- "A"
  } else if ((aggregate[r,"AddRE"] == FALSE) && (aggregate[r,"MulREdim"] > 0)) {
    aggregate[r,"RndEffects"] <- paste("M", aggregate[r,"MulREdim"], sep="_")
  } else {
    aggregate[r,"RndEffects"] <- paste("AM", aggregate[r,"MulREdim"], sep="_")
  }
  aggregate[r,"Var"] <- substr(aggregate[r,"Variable"],1,2)
}

aggregate[,"coverage"] <- aggregate[,"n_covered"] / aggregate[,"n_runs"]
aggregate[,"estimage_sd"] <- sqrt(aggregate[,"estimate_var"])

###############################################################################

pdf(outfile, width=10.5, height=8.5)

# Coverage for data with no hidden covariates: by variable and random effects
ggplot(
  subset(aggregate, (Z_additive == FALSE)&(Z_multiplicative==FALSE)),
  aes(x=RndEffects, y=coverage)
  ) +
  geom_boxplot(aes(col=Variable)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="Coverage with No Unobserved Covariates",
       x="Random Effects", y="Coverage")
# Coverage for data with a multiplicative hidden covariate: by variable and random effects
ggplot(
  subset(aggregate, (Z_additive == FALSE)&(Z_multiplicative==TRUE)),
  aes(x=RndEffects, y=coverage)
) +
  geom_boxplot(aes(col=Variable)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="Coverage with a 2-dimensional Multiplicative Unobserved Covariate",
       x="Random Effects", y="Coverage")


# Bias for data with no hidden covariates: by variable and random effects
ggplot(
  subset(aggregate, (Z_additive == FALSE)&(Z_multiplicative==FALSE)),
  aes(x=RndEffects, y=estimate_mean)
  ) +
  geom_boxplot(aes(col=Variable)) +
  geom_abline(linetype=2, slope=0, intercept=1) +
  geom_abline(linetype=2, slope=0, intercept=-1) +
  labs(title="Bias with No Unobserved Covariates",
       x="Random Effects", y="Mean of Estimates")
# Bias for data with a multiplicative hidden covariate: by variable and random effects
ggplot(
  subset(aggregate, (Z_additive == FALSE)&(Z_multiplicative==TRUE)),
  aes(x=RndEffects, y=estimate_mean)
  ) +
  geom_boxplot(aes(col=Variable)) +
  geom_abline(linetype=2, slope=0, intercept=1) +
  geom_abline(linetype=2, slope=0, intercept=-1) +
  labs(title="Bias with a 2-dimensional Multiplicative Unobserved Covariate",
       x="Random Effects", y="Mean of Estimates")

dev.off()

