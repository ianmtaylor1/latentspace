# Analyze the results of the first Simulation using the custom version of the code
library(RColorBrewer)
library(ggplot2)

infile <- "data/JSM Results n=40.csv"
outfile <- "plots.pdf"

mypalette <- c(
  brewer.pal(10,"Paired")[c(1,3,5,9)]
)

palette(mypalette)


results <- read.csv(infile, header=TRUE, stringsAsFactors=FALSE)

results[,"Covered"] <- 1 * ((results[,"TrueValue"] >= results[,"CI_low"]) & (results[,"TrueValue"] <= results[,"CI_high"]))
results[,"CI_width"] <- results[,"CI_high"] - results[,"CI_low"]

aggregate <- unique(results[,c("Run","Xseed","Z_additive","Z_multiplicative","AddRE","Variable","TrueValue","prior","projected")])
aggregate[,"Prior"] <- "ERR"

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
    aggregate[r,"Prior"] <- ifelse(aggregate[r,"prior"]=="IG","Inverse Wishart",ifelse(aggregate[r,"prior"]=="HC","Half Cauchy","ERR"))
  }
  aggregate[r,"Var"] <- ifelse(aggregate[r,"Variable"] == "intercept", "beta[0]", paste("beta[",substr(aggregate[r,"Variable"],2,2),"]",sep=""))
}

aggregate[,"coverage"] <- aggregate[,"n_covered"] / aggregate[,"n_runs"]
aggregate[,"estimage_sd"] <- sqrt(aggregate[,"estimate_var"])
aggregate$RndEffects <- factor(aggregate$RndEffects, levels=c("None","A_IG","A_HC"))
aggregate$projected <- factor(aggregate$projected, levels=c(FALSE, TRUE))
aggregate$Prior <- factor(aggregate$Prior, levels=c("Inverse Wishart", "Half Cauchy"))



###############################################################################

pdf(outfile, width=9, height=6)


# Coverage for when there is no nodal variation using random effects with default prior
ggplot(
  subset(aggregate, (Z_additive == FALSE)&(RndEffects=="A_IG")&(projected==FALSE)),
  aes(x=Var, y=coverage)
  ) +
  geom_boxplot(aes(col=Prior)) + guides(color=FALSE) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(x="", y="Coverage") +
  theme_bw() +
  theme(text=element_text(size=24), legend.position="bottom") +
  scale_color_manual(values=brewer.pal(10,"Paired")[c(2,6)]) +
  scale_x_discrete(labels=function(x) {parse(text=x)})


dev.off()

