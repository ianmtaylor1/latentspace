# Analyze the results of the first Simulation using the custom version of the code
library(RColorBrewer)
library(ggplot2)
library(here)

jsminfile <- here("jsmposter","data","JSM Results n=40.csv")
outdir <- here("manuscript", "plots")

mypalette <- c(
  brewer.pal(10,"Paired")[c(1,3,5,9)]
)

palette(mypalette)


results <- read.csv(jsminfile, header=TRUE, stringsAsFactors=FALSE)

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

# # Coverage for when there is no nodal variation - priors
# ggplot(
#   subset(aggregate, (Z_additive == FALSE)&(RndEffects!="None")&(projected==FALSE)),
#   aes(x=Var, y=coverage)
# ) +
#   geom_boxplot(aes(col=Prior)) + #guides(color=FALSE) +
#   coord_cartesian(ylim=c(0.5,1)) +
#   geom_abline(linetype=2, slope=0, intercept=.9) +
#   labs(x="", y="Coverage") +
#   theme_bw() +
#   theme(text=element_text(size=30), legend.position=c(0.8, 0.1)) +
#   scale_color_manual(values=brewer.pal(10,"Paired")[c(2,6)]) +
#   scale_x_discrete(labels=function(x) {parse(text=x)})
# ggsave(file.path(outdir, "coverage-novariation.png"), width=9, height=9)
# 
# # Coverage for when there is nodal variation - priors
# ggplot(
#   subset(aggregate, (Z_additive == TRUE)&(RndEffects!="None")&(projected==FALSE)),
#   aes(x=Var, y=coverage)
# ) +
#   geom_boxplot(aes(col=Prior)) + #guides(color=FALSE) +
#   coord_cartesian(ylim=c(0.5,1)) +
#   geom_abline(linetype=2, slope=0, intercept=.9) +
#   labs(x="", y="Coverage") +
#   theme_bw() +
#   theme(text=element_text(size=30), legend.position=c(0.8, 0.1)) +
#   scale_color_manual(values=brewer.pal(10,"Paired")[c(2,6)]) +
#   scale_x_discrete(labels=function(x) {parse(text=x)})
# ggsave(file.path(outdir, "coverage-variation.png"), width=9, height=9)

# Coverage with and without nodal variation - priors
ggplot(
  subset(aggregate, (RndEffects!="None")&(projected==FALSE)),
  aes(x=Var, y=coverage)
) +
  geom_boxplot(aes(col=Prior)) + #guides(color=FALSE) +
  facet_wrap("Z_additive", labeller=labeller(Z_additive=function(x) paste0(ifelse(x, "(B) ", "(A) No "), "Unobserved Nodal Variation"))) +
  coord_cartesian(ylim=c(0.5,1)) +
  geom_abline(linetype=2, slope=0, intercept=.9) +
  labs(x="", y="Coverage", col="Variance Prior") +
  theme_bw() +
  theme(text=element_text(size=36), legend.position=c(0.8, 0.15)) +
  scale_color_manual(values=brewer.pal(10,"Paired")[c(2,6)]) +
  scale_x_discrete(labels=function(x) {parse(text=x)})
ggsave(file.path(outdir, "coverage-with-and-without-variation.png"), width=18, height=9)


