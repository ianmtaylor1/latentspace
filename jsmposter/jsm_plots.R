# Analyze the results of the first Simulation using the custom version of the code
library(RColorBrewer)
library(ggplot2)

setwd("/home/ian/Documents/Git/latentspace/jsmposter")

infile <- "data/JSM Results n=20.csv"
outfile <- "JSM Results n=20.pdf"

mypalette <- c(
  brewer.pal(10,"Paired")[c(1,3,5,9)]
)

palette(mypalette)


results <- read.csv(infile, header=TRUE, stringsAsFactors=FALSE)

results[,"Covered"] <- 1 * ((results[,"TrueValue"] >= results[,"CI_low"]) & (results[,"TrueValue"] <= results[,"CI_high"]))
results[,"CI_width"] <- results[,"CI_high"] - results[,"CI_low"]

aggregate <- unique(results[,c("Run","Xseed","Z_additive","Z_multiplicative","AddRE","Variable","TrueValue","prior")])

for (r in 1:(dim(aggregate)[1])) {
  subresults <- results[(results[,"Run"] == aggregate[r,"Run"]) # Xseed is based on Run
                        & (results[,"Z_additive"] == aggregate[r,"Z_additive"])
                        & (results[,"Z_multiplicative"] == aggregate[r,"Z_multiplicative"])
                        & (results[,"AddRE"] == aggregate[r,"AddRE"])
                        & (results[,"Variable"] == aggregate[r,"Variable"]) # TrueValue is based on Variable
                        & (results[,"prior"] == aggregate[r,"prior"]),]
  aggregate[r,"n_runs"] <- dim(subresults)[1]
  aggregate[r,"n_covered"] <- sum(subresults[,"Covered"])
  aggregate[r,"estimate_mean"] <- mean(subresults[,"Estimate"])
  aggregate[r,"estimate_var"] <- var(subresults[,"Estimate"])
  aggregate[r,"mean_CI_width"] <- mean(subresults[,"CI_width"])
  if((aggregate[r,"AddRE"] == FALSE)) {
    aggregate[r,"RndEffects"] <- "None"
  } else {
    aggregate[r,"RndEffects"] <- paste("A",aggregate[r,"prior"],sep="_")
  }
  aggregate[r,"Var"] <- paste("\\beta_",substr(aggregate[r,"Variable"],2,2),sep="")
}

aggregate[,"coverage"] <- aggregate[,"n_covered"] / aggregate[,"n_runs"]
aggregate[,"estimage_sd"] <- sqrt(aggregate[,"estimate_var"])

###############################################################################

pdf(outfile, width=9, height=6)

# Coverage for data with no hidden covariates and inverse-gammma prior
ggplot(
  subset(aggregate, (Z_additive == FALSE)&(Z_multiplicative==FALSE)&(RndEffects=="A_IG")),
  aes(x=Variable, y=coverage)
  ) +
  geom_boxplot(aes(col=Variable)) + guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="No Nodal Variation, Inverse-Gamma",
       x="", y="Coverage") +
  theme(text=element_text(size=24))

# Coverage for data with no hidden covariates and half-cauchy prior
ggplot(
  subset(aggregate, (Z_additive == FALSE)&(Z_multiplicative==FALSE)&(RndEffects=="A_HC")),
  aes(x=Variable, y=coverage)
) +
  geom_boxplot(aes(col=Variable)) + guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="No Nodal Variation, Half-Cauchy",
       x="", y="Coverage") +
  theme(text=element_text(size=24))

# Coverage for data with hidden covariates and inverse-gammma prior
ggplot(
  subset(aggregate, (Z_additive == TRUE)&(Z_multiplicative==FALSE)&(RndEffects=="A_IG")),
  aes(x=Variable, y=coverage)
) +
  geom_boxplot(aes(col=Variable)) + guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="Nodal Variation, Inverse-Gamma",
       x="", y="Coverage") +
  theme(text=element_text(size=24))

# Coverage for data with hidden covariates and half-cauchy prior
ggplot(
  subset(aggregate, (Z_additive == TRUE)&(Z_multiplicative==FALSE)&(RndEffects=="A_HC")),
  aes(x=Variable, y=coverage)
) +
  geom_boxplot(aes(col=Variable)) + guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="Nodal Variation, Half-Cauchy",
       x="", y="Coverage") +
  theme(text=element_text(size=24))

#######################################################33

# Coverage for when there is no nodal variation
ggplot(
  subset(aggregate, (Z_additive == FALSE)&(RndEffects!="None")),
  aes(x=Var, y=coverage)
  ) +
  geom_boxplot(aes(col=RndEffects)) + #guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="No Nodal Variation",
       x="", y="Coverage") +
  theme(text=element_text(size=24))


# Coverage for when there is nodal variation
ggplot(
  subset(aggregate, (Z_additive == TRUE)&(RndEffects!="None")),
  aes(x=Var, y=coverage)
  ) +
  geom_boxplot(aes(col=RndEffects)) + #guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="Nodal Variation",
       x="", y="Coverage") +
  theme(text=element_text(size=24))



dev.off()

