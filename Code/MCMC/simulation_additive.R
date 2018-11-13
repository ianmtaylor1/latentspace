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

##################################################################################################################

library(truncnorm)

#########################################################################
###### HELPER FUNCTIONS #################################################
#########################################################################

# Returns a boolean array to use as a filter to remove NA's from a 
# columnized matrix with NA's on the diagonal (i.e if Y is an n-by-n
# matrix with NA's on the diagonal, we want to remove the NA's from
# c(Y).)
nafilter <- function(n) {
  filter <- rep(TRUE, n^2)
  for (i in 1:n) {
    filter[n * (i - 1) + i] <- FALSE
  }
  return(filter)
}


#########################################################################
###### FULL CONDITONAL UPDATE FUNCTIONS #################################
#########################################################################

update_beta_fc <- function(Y, Z, X, sigma.squared, a, b, sigma.ab) {
  # Hyperparameters
  n <- dim(Z)[1]
  p <- dim(X)[3]
  Xmat <- apply(X, 3, c)[nafilter(n),]
  XtX <- t(Xmat) %*% Xmat
  Sigma_beta <- solve(XtX) * (n ^ 2)
  # Base matrix used in both mean and covar of beta
  V <- solve(XtX + solve(Sigma_beta / sigma.squared))
  # Y without all other effects
  zprime <- Z - (a %*% t(rep(1, n))) - (rep(1, n) %*% t(b))
  zcol <- c(zprime)[nafilter(n)]
  # The mean of beta's full conditional distribution
  fcmean <- V %*% t(Xmat) %*% zcol
  # Generate new beta
  beta <- c(rmvnorm(n=1, mu=fcmean, Sigma=sigma.squared * V))
  return(beta)
}

update_sigma.squared_fc <- function(Y, Z, X, beta, a, b, sigma.ab) {
  # Hyperparameters
  alpha_sigma <- 0.5
  beta_sigma <- 0.5
  # Find n
  n <- dim(Z)[1]
  # Y without all other effects
  zprime <- Z - (a %*% t(rep(1, n))) - (rep(1, n) %*% t(b)) - Xbeta(X, beta)
  zcol <- c(zprime)[nafilter(n)]
  # Parameters of fc gamma distribution of sigma.squared inverse
  shape <- ((n * (n - 1) / 2) + alpha_sigma)
  rate <- (t(zcol) %*% zcol) / 2 + 1/beta_sigma
  # Generate new sigma.squared
  sigma.squared <- 1/rgamma(n=1, shape=shape, rate=rate)
  return(sigma.squared)
}

update_a_fc <- function(Y, Z, X, beta, sigma.squared, b, sigma.ab) {
  # Calculate rowmeans of Y, without other effects
  zprime <- Z - (rep(1, n) %*% t(b)) - Xbeta(X, beta)
  n <- dim(Z)[1]
  zrowmeans <- rowSums(zprime, na.rm=TRUE) / (n - 1)
  # Calculate variance and meanvector of 
  var.a <- 1 / (((n - 1) / sigma.squared) + (sigma.ab[2,2] / det(sigma.ab)))
  mean.a <- var.a * (((n-1) * zrowmeans / sigma.squared) + (sigma.ab[1,2] * b / det(sigma.ab)))
  # Calculate new a's
  a <- rnorm(n=n, mean=mean.a, sd=sqrt(var.a))
  return(a)
}

update_b_fc <- function(Y, Z, X, beta, sigma.squared, a, sigma.ab) {
  # Calculate colmeans of Y, without other effects
  zprime <- Z - (a %*% t(rep(1, n))) - Xbeta(X, beta)
  n <- dim(Z)[1]
  zcolmeans <- colSums(zprime, na.rm=TRUE) / (n - 1)
  # Calculate variance and meanvector of 
  var.b <- 1 / (((n - 1) / sigma.squared) + (sigma.ab[1,1] / det(sigma.ab)))
  mean.b <- var.b * (((n-1) * zcolmeans / sigma.squared) + (sigma.ab[1,2] * a / det(sigma.ab)))
  # Calculate new b's
  b <- rnorm(n=n, mean=mean.b, sd=sqrt(var.b))
  return(b)
}

update_sigma.ab_fc <- function(Y, Z, X, beta, sigma.squared, a, b) {
  # Hyperparameters
  V_sigma.ab <- diag(c(1,1))
  n_sigma.ab <- 3
  # Find n
  n <- dim(Y)[1]
  # Create matrix from a and b for use in new parameters
  abmat <- matrix(c(
    sum(a^2),
    sum(a*b),
    sum(a*b),
    sum(b^2)
  ), ncol=2)
  # New parameters
  Snew <- solve(solve(V_sigma.ab) + abmat)
  nunew <- n_sigma.ab + n
  # Calculate the new value and return
  return(solve(rwish(S0=Snew, nu=nunew)))
}

update_Z_bin_fc <- function(Y, X, beta, sigma.squared, a, b, sigma.ab) {
  # In binary (probit) case, the full conditional of z is a truncated normal
  n <- dim(Y)[1]
  Z <- matrix(rep(NA,n*n), ncol=n)
  Zmean <- c(Xbeta(X, beta) + (a %*% t(rep(1, n))) + (rep(1, n) %*% t(b)))
  highbound <- replace(x=rep(0, n^2), (c(Y) == 1), Inf)
  lowbound <- replace(rep(0, n^2), (c(Y) == 0), -Inf)
  Zcol <- rtruncnorm(n=n^2, a=lowbound, b=highbound, mean=Zmean, sd=sqrt(sigma.squared))
  Z <- matrix(Zcol, ncol=n)
  diag(Z) <- NA
  return(Z)
}
#########################################################################
###### MAIN MCMC FUNCTION ###############################################
#########################################################################

imt_ame <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, intercept=TRUE,
                    model="nrm", rvar=TRUE, cvar=TRUE,
                    burn=500, n.iter=10000, save.every=25, print=TRUE) {
  # Check for appropriate model
  stopifnot(model %in% c("nrm", "bin"))
  
  # Set up data as appropriate
  stopifnot(length(dim(Y)) == 2)
  stopifnot(dim(Y)[1] == dim(Y)[2])
  n <- dim(Y)[1]
  X <- design_array(Xrow, Xcol, Xdyad, intercept, n)
  p <- dim(X)[3]
  diag(Y) <- NA
  
  # Initial values
  a <- rep(0, n)
  b <- rep(0, n)
  beta <- rep(0, p)
  sigma.squared <- 1
  sigma.ab <- diag(2)
  if (model == "nrm") {
    Z <- Y
  } else if (model == "bin") {
    Z <- matrix(rep(0,n*n), ncol=n)
    diag(Z) <- NA
  }
  
  # Lists to save the samples in temporarily
  samples.a <- list()
  samples.b <- list()
  samples.beta <- list()
  samples.sigma.squared <- list()
  samples.sigma.ab <- list()
  
  # Iterate through the gibbs sampler
  for (i in 1:(n.iter + burn)) {
    # Update latent normal variables (if probit model)
    if (model == "nrm") {
      Z <- Y
    } else if (model == "bin") {
      Z <- update_Z_bin_fc(Y=Y, X=X, beta=beta, sigma.squared=sigma.squared, a=a, b=b, sigma.ab=sigma.ab)
    }
    # Update sigma.squared
    if (model == "nrm") {
      sigma.squared <- update_sigma.squared_fc(Y=Y, Z=Z, X=X, beta=beta, a=a, b=b, sigma.ab=sigma.ab)
    } else if (model == "bin") {
      sigma.squared <- 1
    }
    # Update Sigma.ab
    if (rvar || cvar) {
      sigma.ab <- update_sigma.ab_fc(Y=Y, Z=Z, X=X, beta=beta, sigma.squared=sigma.squared, a=a, b=b)
    }
    # Update beta
    beta <- update_beta_fc(Y=Y, Z=Z, X=X, sigma.squared=sigma.squared, a=a, b=b, sigma.ab=sigma.ab)
    # Update a,b
    if (rvar) {
      a <- update_a_fc(Y=Y, Z=Z, X=X, beta=beta, sigma.squared=sigma.squared, b=b, sigma.ab=sigma.ab)
    }
    if (cvar) {
      b <- update_b_fc(Y=Y, Z=Z, X=X, beta=beta, sigma.squared=sigma.squared, a=a, sigma.ab=sigma.ab)
    }
    
    # Save if necessary
    if ((i > burn) && ((i-burn) %% save.every == 0)) {
      if (print) {
        cat("Saving...",i,"\n")
      }
      samples.a <- append(samples.a, list(a))
      samples.b <- append(samples.b, list(b))
      samples.beta <- append(samples.beta, list(beta))
      samples.sigma.squared <- append(samples.sigma.squared, list(sigma.squared))
      samples.sigma.ab <- append(samples.sigma.ab, list(sigma.ab))
    }
  }
  
  # Create form in which to return results
  results.beta <- matrix(unlist(samples.beta), ncol=p, byrow=TRUE)
  colnames(results.beta) <- dimnames(X)[[3]]
  results.a <- matrix(unlist(samples.a), ncol=n, byrow=TRUE)
  results.b <- matrix(unlist(samples.b), ncol=n, byrow=TRUE)
  vcdata <- c()
  for (i in 1:length(samples.sigma.ab)) {
    vcdata <- c(vcdata, 
                samples.sigma.ab[[i]][1,1], #va
                samples.sigma.ab[[i]][2,2], #vb
                samples.sigma.ab[[i]][1,2], #cab
                0, #rho
                samples.sigma.squared[[i]] #ve
    )
  }
  results.vc <- matrix(vcdata, ncol=5, byrow=TRUE)
  colnames(results.vc) <- c("va","vb","cab","rho","ve")
  
  # Return Samples
  return(
    list(
      BETA=results.beta,
      A=results.a,
      B=results.b,
      VC=results.vc
    )
  )
}

##################################################################################################################



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
# ADDITIVE ROW AND COLUMN EFFECTS
true.a <- rnorm(n=n)
true.b <- rnorm(n=n)
Z <- outer(true.a, true.b, "+")
diag(Z) <- NA

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
  arglists <- list(list(Y=Y_binary_noZ, Z=FALSE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter),
                   list(Y=Y_binary, Z=TRUE, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter),
                   list(Y=Y_binary_noZ, Z=FALSE, Xdyad=Xd, rvar=FALSE, cvar=FALSE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter),
                   list(Y=Y_binary, Z=TRUE, Xdyad=Xd, rvar=FALSE, cvar=FALSE, dcor=FALSE, R=0, plot=FALSE, print=FALSE, model="bin", gof=FALSE, nscan=mcmc.iter)
  )
  
  # Run the model once for each set of arguments
  for (args in arglists) {
    fit <- do.call(imt_ame, args[c("Y","Xdyad","model","rvar","cvar","print")])
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


