# Analyze the results of the first Simulation using the custom version of the code
library(RColorBrewer)

mypalette <- c(
  brewer.pal(10,"Paired")[c(1,3,5,9)]
)

palette(mypalette)


results <- read.csv("Results/Global_Horseshoe_prior_results.csv", header=TRUE, stringsAsFactors=FALSE)

results[,"Covered"] <- 1 * ((results[,"TrueValue"] >= results[,"CI_low"]) & (results[,"TrueValue"] <= results[,"CI_high"]))
results[,"CI_width"] <- results[,"CI_high"] - results[,"CI_low"]

aggregate <- unique(results[,c("Run","Xseed","Z_additive","Z_multiplicative","AddRE","MulRE","Variable","TrueValue","prior")])

for (r in 1:(dim(aggregate)[1])) {
  subresults <- results[(results[,"Run"] == aggregate[r,"Run"]) # Xseed is based on Run
                        & (results[,"Z_additive"] == aggregate[r,"Z_additive"])
                        & (results[,"Z_multiplicative"] == aggregate[r,"Z_multiplicative"])
                        & (results[,"AddRE"] == aggregate[r,"AddRE"])
                        & (results[,"MulRE"] == aggregate[r,"MulRE"])
                        & (results[,"Variable"] == aggregate[r,"Variable"]),] # TrueValue is based on Variable
  aggregate[r,"n_runs"] <- dim(subresults)[1]
  aggregate[r,"n_covered"] <- sum(subresults[,"Covered"])
  aggregate[r,"estimate_mean"] <- mean(subresults[,"Estimate"])
  aggregate[r,"estimate_var"] <- var(subresults[,"Estimate"])
  aggregate[r,"mean_CI_width"] <- mean(subresults[,"CI_width"])
  if((aggregate[r,"AddRE"] == FALSE) && (aggregate[r,"MulRE"] == FALSE)) {
    aggregate[r,"RndEffects"] <- ""
  } else if ((aggregate[r,"AddRE"] == TRUE) && (aggregate[r,"MulRE"] == FALSE)) {
    aggregate[r,"RndEffects"] <- paste("A",aggregate[r,"prior"],sep="_")
  } else if ((aggregate[r,"AddRE"] == FALSE) && (aggregate[r,"MulRE"] == TRUE)) {
    aggregate[r,"RndEffects"] <- "M"
  } else {
    aggregate[r,"RndEffects"] <- "MA"
  }
  aggregate[r,"Var"] <- substr(aggregate[r,"Variable"],1,2)
}

aggregate[,"coverage"] <- aggregate[,"n_covered"] / aggregate[,"n_runs"]
aggregate[,"estimage_sd"] <- sqrt(aggregate[,"estimate_var"])


# Begin creating plots
pdf('rplot.pdf', width=10.5, height=8.5)

# Plot coverage of confidence intervals
boxplot(coverage~Var*RndEffects, data=aggregate[(aggregate[,"Z_additive"]==TRUE),],
        col=c(1,2,3,4),
        main="Data with Additive Unobserved Covariates", ylab="Coverage")
abline(h=0.9,col="Black",lty=2) # HARD CODED
boxplot(coverage~Var*RndEffects, data=aggregate[(aggregate[,"Z_multiplicative"]==TRUE),],
        col=c(1,2,3,4),
        main="Data With Multiplicative Unobserved Covariate", ylab="Coverage")
abline(h=0.9,col="Black",lty=2) # HARD CODED
boxplot(coverage~Var*RndEffects, data=aggregate[(aggregate[,"Z_multiplicative"]==FALSE)&(aggregate[,"Z_additive"]==FALSE),],
        col=c(1,2,3,4),
        main="Data Without Unobserved Covariates", ylab="Coverage")
abline(h=0.9,col="Black",lty=2) # HARD CODED

# Plot means of estimates over all Reps
boxplot(estimate_mean~Var*RndEffects, data=aggregate[(aggregate[,"Z_additive"]==TRUE),],
        col=c(1,2,3,4),
        main="Data with Additive Unobserved Covariates", ylab="Mean of Estimated Coefficients")
abline(h=1,col="Black",lty=2) # HARD CODED
abline(h=-1,col="Black",lty=2) # HARD CODED
boxplot(estimate_mean~Var*RndEffects, data=aggregate[(aggregate[,"Z_multiplicative"]==TRUE),],
        col=c(1,2,3,4),
        main="Data WIith Multiplicative Unobserved Covariate", ylab="Mean of Estimated Coefficients")
abline(h=1,col="Black",lty=2) # HARD CODED
abline(h=-1,col="Black",lty=2) # HARD CODED
boxplot(estimate_mean~Var*RndEffects, data=aggregate[(aggregate[,"Z_additive"]==FALSE)&(aggregate[,"Z_multiplicative"]==FALSE),],
        col=c(1,2,3,4),
        main="Data without Unobserved Covariates", ylab="Mean of Estimated Coefficients")
abline(h=1,col="Black",lty=2) # HARD CODED
abline(h=-1,col="Black",lty=2) # HARD CODED

# Plot variances of estimates over all Reps
boxplot(estimate_var~Var*RndEffects, data=aggregate[(aggregate[,"Z_additive"]==TRUE),],
        col=c(1,2,3,4),
        main="Data with Additive Unobserved Covariates", ylab="Variance of Estimated Coefficients")
boxplot(estimate_var~Var*RndEffects, data=aggregate[(aggregate[,"Z_multiplicative"]==TRUE),],
        col=c(1,2,3,4),
        main="Data With Multiplicative Unobserved Covariate", ylab="Variance of Estimated Coefficients")
boxplot(estimate_var~Var*RndEffects, data=aggregate[(aggregate[,"Z_additive"]==FALSE)&(aggregate[,"Z_multiplicative"]==FALSE),],
        col=c(1,2,3,4),
        main="Data without Unobserved Covariate", ylab="Variance of Estimated Coefficients")

# Plot average CI width over all Reps
boxplot(mean_CI_width~Var*RndEffects, data=aggregate[(aggregate[,"Z_additive"]==TRUE),],
        col=c(1,2,3,4),
        main="Data with Additive Unobserved Covariates", ylab="Average width of Confidence Interval")
boxplot(mean_CI_width~Var*RndEffects, data=aggregate[(aggregate[,"Z_multiplicative"]==TRUE),],
        col=c(1,2,3,4),
        main="Data With Multiplicative Unobserved Covariate", ylab="Average width of Confidence Interval")
boxplot(mean_CI_width~Var*RndEffects, data=aggregate[(aggregate[,"Z_additive"]==FALSE)&(aggregate[,"Z_multiplicative"]==FALSE),],
        col=c(1,2,3,4),
        main="Data without Unobserved Covariates", ylab="Average width of Confidence Interval")

dev.off()

