# Analyze the results of the first Simulation using the custom version of the code
library(RColorBrewer)
library(ggplot2)

setwd("/home/ian/Documents/Git/latentspace/jsmposter")

infile <- "data/JSM Results n=80.csv"
outfile <- "JSM Results n=80.pdf"

mypalette <- c(
  brewer.pal(10,"Paired")[c(1,3,5,9)]
)

palette(mypalette)


results <- read.csv(infile, header=TRUE, stringsAsFactors=FALSE)

results[,"Covered"] <- 1 * ((results[,"TrueValue"] >= results[,"CI_low"]) & (results[,"TrueValue"] <= results[,"CI_high"]))
results[,"CI_width"] <- results[,"CI_high"] - results[,"CI_low"]

aggregate <- unique(results[,c("Run","Xseed","Z_additive","Z_multiplicative","AddRE","Variable","TrueValue","prior","projected")])

for (r in 1:(dim(aggregate)[1])) {
  subresults <- results[(results[,"Run"] == aggregate[r,"Run"]) # Xseed is based on Run
                        & (results[,"Z_additive"] == aggregate[r,"Z_additive"])
                        & (results[,"Z_multiplicative"] == aggregate[r,"Z_multiplicative"])
                        & (results[,"AddRE"] == aggregate[r,"AddRE"])
                        & (results[,"Variable"] == aggregate[r,"Variable"]) # TrueValue is based on Variable
                        & (results[,"prior"] == aggregate[r,"prior"])
			& (results[,"projected"] == aggregate[r,"projected"]),]
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
  aggregate[r,"Var"] <- ifelse(aggregate[r,"Variable"] == "intercept", "\\beta_0", paste("\\beta_",substr(aggregate[r,"Variable"],2,2),sep=""))
}

aggregate[,"coverage"] <- aggregate[,"n_covered"] / aggregate[,"n_runs"]
aggregate[,"estimage_sd"] <- sqrt(aggregate[,"estimate_var"])
aggregate$RndEffects <- factor(aggregate$RndEffects, levels=c("None","A_IG","A_HC"))
aggregate$projected <- factor(aggregate$projected, levels=c(FALSE, TRUE))



###############################################################################

pdf(outfile, width=9, height=6)


# Coverage for when there is no nodal variation - priors
ggplot(
  subset(aggregate, (Z_additive == FALSE)&(RndEffects!="None")&(projected==FALSE)),
  aes(x=Var, y=coverage)
  ) +
  geom_boxplot(aes(col=RndEffects)) + guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="Prior Comparison: No Nodal Variation",
       x="", y="Coverage") +
  theme(text=element_text(size=24)) +
  scale_color_manual(values=brewer.pal(10,"Paired")[c(1,3)])


# Coverage for when there is nodal variation - priors
ggplot(
  subset(aggregate, (Z_additive == TRUE)&(RndEffects!="None")&(projected==FALSE)),
  aes(x=Var, y=coverage)
  ) +
  geom_boxplot(aes(col=RndEffects)) + guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="Prior Comparison: Nodal Variation",
       x="", y="Coverage") +
  theme(text=element_text(size=24)) +
  scale_color_manual(values=brewer.pal(10,"Paired")[c(1,3)])

# Coverage for when there is no nodal variation - projections
ggplot(
  subset(aggregate, (Z_additive == FALSE)&(RndEffects!="None")&(RndEffects=="A_IG")),
  aes(x=Var, y=coverage)
  ) +
  geom_boxplot(aes(col=projected)) + guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="Projection Comparison: No Nodal Variation",
       x="", y="Coverage") +
  theme(text=element_text(size=24)) +
  scale_color_manual(values=brewer.pal(10,"Paired")[c(1,5)])


# Coverage for when there is nodal variation - projections
ggplot(
  subset(aggregate, (Z_additive == TRUE)&(RndEffects!="None")&(RndEffects=="A_IG")),
  aes(x=Var, y=coverage)
  ) +
  geom_boxplot(aes(col=projected)) + guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(title="Projection Comparison: Nodal Variation",
       x="", y="Coverage") +
  theme(text=element_text(size=24)) +
  scale_color_manual(values=brewer.pal(10,"Paired")[c(1,5)])



dev.off()

