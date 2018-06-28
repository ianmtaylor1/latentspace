library(amen)
library(doParallel)
library(batch)

# Set default values of parameters passed via command line
jobname <- "" # How to name the output file
n <- 20       # Number of nodes in the network
run <- 1      # Which run to do
reps <- 200   # The number of reps in this run
startrep <- 1 # The rep number to start on
workers <- 1  # Number of parallel workers (cores) to use
# Overwrite them with whatever the user provides, if anything
parseCommandArgs()

reps <- startrep:(reps + startrep - 1) # Number of times to generate Y values for each set of covariates

mcmc.iter <- 10000 # Number of iterations of the MCMC algorithm to do (after "burn") (default: 10000)

intercept <- -1
beta <- c(1,1,1) # True coefficients (not including intercept)
gamma <- 1       # Coefficient for unobserved data

conf <- .90      # Confidence level of intervals to measure

rndseedbase <- 10000000

filename <- paste("latent_var_simulation_results_",jobname,"_",format(Sys.time(),"%Y-%m-%d_%H-%M"),".csv",sep="")


###########################################################
#### Start the run
###########################################################

cat("Run",run,"-",date(),"\n")
  
# ame() resets random seed, so seed each loop with new number
xseed <- rndseedbase * run
set.seed(xseed)

#### Generate design matrix X for dyad covariates
# Start empty
Xd <- array(0, dim=c(n,n,length(beta)) )
# Fill first covariate with binary indicator of shared class
groups <- sample(0:1, n, TRUE) * 2 - 1
if (all(groups == 1) || all(groups == -1)) {
  k <- sample(1:n, 1)[1]
  groups[k] <- -groups[k]
}
Xd[,,1] <- outer(groups, groups, "*")
# Fill second covariate with "distances" between normally distributed coordinates
locations <- rnorm(n=n, mean=0, sd=1)
dist <- function(x,y) {
  abs(x - y)
}
Xd[,,2] <- outer(locations, locations, "dist")
# Fill third covariate with all random standard normal draws
Xd[,,3] <- matrix(rnorm(n=n*n, mean=0, sd=1), nrow=n, ncol=n)
# Set all diagonals to NA
for (i in 1:n) {
  for (k in 1:length(beta)) {
    Xd[i,i,k] <- NA
  }
}

dimnames(Xd)[[3]] <- c("X1","X2","X3")

#### Generate unobserved covariate
# Similar to the shared class but with three classes (b/c why not)
uo_groups <- sample(0:1, n, TRUE) * 2 - 1
if (all(uo_groups == -1) || all(uo_groups == 1)) {
  k <- sample(1:n, 1)[1]
  uo_groups[k] <- -uo_groups[k]
X}
Z <- outer(uo_groups, uo_groups, "*")
for (i in 1:n) {
  Z[i,i] <- NA
}

########################################################
#### Loop through each rep in parallel
########################################################

registerDoParallel(workers, cores=workers)

results <- foreach(rep=reps, .combine="rbind") %dopar% {
  cat("Run",run,"Rep",rep,"-",date(),"\n")
  
  library(amen)
  
  represults <- data.frame(Run=integer(),
                           Rep=integer(),
                           Xseed=integer(),
                           Yseed=integer(),
                           Z=logical(),
                           AddRE=logical(),
                           MulRE=logical(),
                           Variable=character(),
                           TrueValue=double(),
                           Estimate=double(),
                           Variance=double(),
                           Confidence=double(),
                           CI_low=double(),
                           CI_high=double(),
                           stringsAsFactors=FALSE)
  
  yseed <- rndseedbase * run + rep
  set.seed(yseed)

  #### Simulate binary data with Probita
  eta <- intercept
  for (k in 1:length(beta)) {
    eta <- eta + beta[k]*Xd[,,k]
  }
  eta <- eta + matrix(rnorm(n = n*n), nrow=n, ncol=n)
  Y_binary_noZ <- 1 * (eta > 0) # Data generated WITHOUT unobserved covariate
  Y_binary <- 1 * (eta + Z * gamma > 0) # Data generated WITH unobserved covariate
  
  #### Fit models and record confidence intervals
  arglists <- list(list(Y=Y_binary_noZ, Z=FALSE, Xdyad=Xd, rvar=FALSE, cvar=FALSE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter),
                   list(Y=Y_binary_noZ, Z=FALSE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter),
                   list(Y=Y_binary_noZ, Z=FALSE, Xdyad=Xd, rvar=FALSE, cvar=FALSE, dcor=FALSE, R=1, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter),
                   list(Y=Y_binary_noZ, Z=FALSE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=1, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter),
                   list(Y=Y_binary, Z=TRUE, Xdyad=Xd, rvar=FALSE, cvar=FALSE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter),
                   list(Y=Y_binary, Z=TRUE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter),
                   list(Y=Y_binary, Z=TRUE, Xdyad=Xd, rvar=FALSE, cvar=FALSE, dcor=FALSE, R=1, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter),
                   list(Y=Y_binary, Z=TRUE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=1, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter)
  )
  
  # Run the model once for each set of arguments
  for (args in arglists) {
    fit <- do.call(ame, args[c("Y","Xdyad","rvar","cvar","dcor","R","plot","print","model","gof","nscan")])
    variables <- list(list(name="intercept", value=intercept, fitname="intercept"),
                      list(name="X1_groups", value=beta[1], fitname="X1.dyad"),
                      list(name="X2_dist", value=beta[2], fitname="X2.dyad"),
                      list(name="X3_norm", value=beta[3], fitname="X3.dyad"))
    for (v in variables) {
      varrow <- data.frame(Run=run,
                           Rep=rep,
                           Xseed=xseed,
                           Yseed=yseed,
                           Z=args$Z,
                           AddRE=args$rvar,
                           MulRE=(args$R > 0),
                           Variable=v$name,
                           TrueValue=v$value,
                           Estimate=mean(fit$BETA[,v$fitname]),
                           Variance=var(fit$BETA[,v$fitname]),
                           Confidence=conf,
                           CI_low=quantile(fit$BETA[,v$fitname],(1-conf)/2),
                           CI_high=quantile(fit$BETA[,v$fitname],1-(1-conf)/2),
                           stringsAsFactors=FALSE
                           )
      represults <- rbind(represults, varrow)
    }
  }
  # Return the results for this rep
  represults
}

stopImplicitCluster()

# Output the results
write.csv(results, filename, row.names=FALSE)


