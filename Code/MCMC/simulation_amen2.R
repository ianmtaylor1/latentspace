library(amen2)
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
beta <- c(1,1,1,1) # True coefficients (not including intercept)
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
# Fill fourth covariate with a random row effect
Xd[,,4] <- matrix(rnorm(n=n, mean=0, sd=1), nrow=n, ncol=n)
# Set all diagonals to NA
for (i in 1:n) {
  for (k in 1:length(beta)) {
    Xd[i,i,k] <- NA
  }
}

dimnames(Xd)[[3]] <- c("X1","X2","X3","X4")

#### Generate unobserved covariate
# ADDITIVE ROW AND COLUMN EFFECTS
true.a <- rnorm(n=n)
true.b <- rnorm(n=n)
Z.additive <- outer(true.a, true.b, "+")
diag(Z.additive) <- NA

#### Generate unobserved covariate
# Similar to the shared class but with three classes (b/c why not)
uo_groups <- sample(0:1, n, TRUE) * 2 - 1
if (all(uo_groups == -1) || all(uo_groups == 1)) {
  k <- sample(1:n, 1)[1]
  uo_groups[k] <- -uo_groups[k]
  X}
Z.multiplicative <- outer(uo_groups, uo_groups, "*")
diag(Z.multiplicative) <- NA



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
                           Z_additive=logical(),
                           Z_multiplicative=logical(),
                           prior=character(),
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
  Y_binary_Zadditive <- 1 * (eta + Z.additive * gamma > 0) # Data generated WITH unobserved covariate
  Y_binary_Zmultiplicative <- 1 * (eta + Z.multiplicative * gamma > 0) # Data generated WITH unobserved covariate
  
  #### Fit models and record confidence intervals
  arglists <- list(
    # No rvar or cvar
    list(Y=Y_binary_noZ, Z_additive=FALSE, Z_multiplicative=FALSE, Xdyad=Xd, rvar=FALSE, cvar=FALSE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="ig"),
    list(Y=Y_binary_Zadditive, Z_additive=TRUE, Z_multiplicative=FALSE, Xdyad=Xd, rvar=FALSE, cvar=FALSE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="ig"),
    list(Y=Y_binary_Zmultiplicative, Z_additive=FALSE, Z_multiplicative=TRUE, Xdyad=Xd, rvar=FALSE, cvar=FALSE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="ig"),
    # ig prior
    list(Y=Y_binary_noZ, Z_additive=FALSE, Z_multiplicative=FALSE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="ig"),
    list(Y=Y_binary_Zadditive, Z_additive=TRUE, Z_multiplicative=FALSE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="ig"),
    list(Y=Y_binary_Zmultiplicative, Z_additive=FALSE, Z_multiplicative=TRUE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="ig"),
    # hs prior
    list(Y=Y_binary_noZ, Z_additive=FALSE, Z_multiplicative=FALSE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="hs"),
    list(Y=Y_binary_Zadditive, Z_additive=TRUE, Z_multiplicative=FALSE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="hs"),
    list(Y=Y_binary_Zmultiplicative, Z_additive=FALSE, Z_multiplicative=TRUE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="hs"),
    # tiny prior
    list(Y=Y_binary_noZ, Z_additive=FALSE, Z_multiplicative=FALSE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="tiny"),
    list(Y=Y_binary_Zadditive, Z_additive=TRUE, Z_multiplicative=FALSE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="tiny"),
    list(Y=Y_binary_Zmultiplicative, Z_additive=FALSE, Z_multiplicative=TRUE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter, prior="tiny")
  )
  
  # Run the model once for each set of arguments
  for (args in arglists) {
    fit <- do.call(ame2, args[c("Y","Xdyad","model","rvar","cvar","print","prior")])
    variables <- list(list(name="intercept", value=intercept, fitname="intercept"),
                      list(name="X1_groups", value=beta[1], fitname="X1.dyad"),
                      list(name="X2_dist", value=beta[2], fitname="X2.dyad"),
                      list(name="X3_norm", value=beta[3], fitname="X3.dyad"),
                      list(name="X4_row", value=beta[4], fitname="X4.dyad"))
    for (v in variables) {
      varrow <- data.frame(Run=run,
                           Rep=rep,
                           Xseed=xseed,
                           Yseed=yseed,
                           Z_additive=args$Z_additive,
                           Z_multiplicative=args$Z_multiplicative,
                           prior=args$prior,
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


