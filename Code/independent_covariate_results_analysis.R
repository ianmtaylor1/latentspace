# Analyze the results of the first Independent Covariate simulation

results <- read.csv("Results/Independent_Covariate_2_nodcor.csv", header=TRUE, stringsAsFactors=FALSE)

results[,"Covered"] <- 1 * ((results[,"TrueValue"] >= results[,"CI_low"]) & (results[,"TrueValue"] <= results[,"CI_high"]))
results[,"CI_width"] <- results[,"CI_high"] - results[,"CI_low"]

results[1:10,"Covered"]

aggregate <- unique(results[,c("Run","Xseed","Z","AddRE","MulRE","Variable","TrueValue")])

for (r in 1:(dim(aggregate)[1])) {
  subresults <- results[(results[,"Run"] == aggregate[r,"Run"]) # Xseed is based on Run
                        & (results[,"Z"] == aggregate[r,"Z"])
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
    aggregate[r,"RndEffects"] <- "A"
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
boxplot(coverage~Var*RndEffects, data=aggregate[(aggregate[,"Z"]==FALSE),],
        col=c("Gray","Red","Green","Blue"),
        main="Data without Unobserved Covariate", ylab="Coverage")
abline(h=0.9,col="Black",lty=2) # HARD CODED
boxplot(coverage~Var*RndEffects, data=aggregate[(aggregate[,"Z"]==TRUE),],
        col=c("Gray","Red","Green","Blue"),
        main="Data WITH Unobserved Covariate", ylab="Coverage")
abline(h=0.9,col="Black",lty=2) # HARD CODED

# Plot means of estimates over all Reps
boxplot(estimate_mean~Var*RndEffects, data=aggregate[(aggregate[,"Z"]==FALSE),],
        col=c("Gray","Red","Green","Blue"),
        main="Data without Unobserved Covariate", ylab="Mean of Estimated Coefficients")
abline(h=1,col="Black",lty=2) # HARD CODED
abline(h=-1,col="Black",lty=2) # HARD CODED
boxplot(estimate_mean~Var*RndEffects, data=aggregate[(aggregate[,"Z"]==TRUE),],
        col=c("Gray","Red","Green","Blue"),
        main="Data WITH Unobserved Covariate", ylab="Mean of Estimated Coefficients")
abline(h=1,col="Black",lty=2) # HARD CODED
abline(h=-1,col="Black",lty=2) # HARD CODED

# Plot variances of estimates over all Reps
boxplot(estimate_var~Var*RndEffects, data=aggregate[(aggregate[,"Z"]==FALSE),],
        col=c("Gray","Red","Green","Blue"),
        main="Data without Unobserved Covariate", ylab="Variance of Estimated Coefficients")
boxplot(estimate_var~Var*RndEffects, data=aggregate[(aggregate[,"Z"]==TRUE),],
        col=c("Gray","Red","Green","Blue"),
        main="Data WITH Unobserved Covariate", ylab="Variance of Estimated Coefficients")

# Plot average CI width over all Reps
boxplot(mean_CI_width~Var*RndEffects, data=aggregate[(aggregate[,"Z"]==FALSE),],
        col=c("Gray","Red","Green","Blue"),
        main="Data without Unobserved Covariate", ylab="Average width of Confidence Interval")
boxplot(mean_CI_width~Var*RndEffects, data=aggregate[(aggregate[,"Z"]==TRUE),],
        col=c("Gray","Red","Green","Blue"),
        main="Data WITH Unobserved Covariate", ylab="Average width of Confidence Interval")


dev.off()

