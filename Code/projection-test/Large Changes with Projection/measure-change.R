
library(amenhs)
library(batch)
library(foreach)

# Takes a matrix X and a vector Y and rescales Y in terms of its orthogonal components
# Until it has the desired canonical correlation with X. Assumes no centering needs to be done
# The returned vector will have the same length/magnitude as the input vector Y
gencancor <- function(X, Y, rho) {
  PX <- X %*% MASS::ginv(t(X) %*% X) %*% t(X)
  Ypar <- c(PX %*% Y)
  Yperp <- Y - Ypar
  Ynew <- (rho * sqrt(sum(Yperp^2)) * Ypar + sqrt(1-rho^2) * sqrt(sum(Ypar^2)) * Yperp)
  return(Ynew * sqrt(sum(Y^2)) / sqrt(sum(Ynew^2)))
}

################################################################################

# Parse command line arguments into the script
# First: set the defaults
run <- 5         # Which run # is this? 1,2,3,...
reps <- 20       # How many times to generate Y and fit?
#pr <- 1          # How many fixed row covariates?
#pc <- 1          # How many fixed column covariates?
noise <- 1       # Noise level (magnitude of unobserved effects relative to cumulative effect of fixed)
#set.corr <- 0.5  # Canonical correlation of unobserved effects with fixed effects
n <- 27          # Number of nodes in networks
outfilebase <- NULL # Base of name of CSV to write to.
# Read the command line. May override variables set in defaults.
parseCommandArgs()

################################################################################

# Parameters to use for each run
pr.all <-   c(1,   5,   5,   1,   5,   5)
pc.all <-   c(1,   1,   5,   1,   1,   5)
corr.all <- c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9)

# Choose based on run
pr <- pr.all[(run-1) %% 6 + 1]
pc <- pc.all[(run-1) %% 6 + 1]
set.corr <- corr.all[(run-1) %% 6 + 1]

# Repeatability
base.seed <- 111

set.seed(base.seed * 1000 + run)

# Parameters for generating data
intercept <- -1
beta.r <- rep(1, pr)
beta.c <- rep(1, pc)

# Generating the fixed things for this run
Xr <- matrix(rnorm(pr*n), ncol=pr)
Xc <- matrix(rnorm(pc*n), ncol=pc)
a.true <- gencancor(cbind(rep(1,n), Xr), rnorm(n), set.corr) * noise * sqrt(sum(beta.r^2))
b.true <- gencancor(cbind(rep(1,n), Xc), rnorm(n), set.corr) * noise * sqrt(sum(beta.c^2))

# MCMC controlling parameters
iter<-100000
burn<-5000
thin<-50


output <- foreach(rep=1:reps, .combine="rbind") %do% {
  cat("Rep",rep,"\n")
  
  set.seed((base.seed * 1000 + run) * 10000 + rep)
  
  # Response variable
  Y <- (intercept + outer(c(Xr %*% beta.r), c(Xc %*% beta.c), "+") + outer(a.true, b.true, "+") + matrix(rnorm(n*n), nrow=n) >= 0)
  
  # Run the MCMC
  res <- ame(Y, Xrow=Xr, Xcol=Xc, family="bin", halfcauchy=TRUE, project=TRUE, dcor=FALSE,
             nscan=iter, burn=burn, odens=thin,
             print=FALSE, plot=FALSE, gof=FALSE)
  
  # Analyze results
  ret <- data.frame(run=run, rep=rep, pr=pr, pc=pc, noise=noise, corr=set.corr, n=n, stringsAsFactors=FALSE)
  # 1. Unprojected estimates and CI's
  ret[["beta_intercept_mean"]] <- mean(res$BETA[,1])
  ret[["beta_intercept_5q"]] <- quantile(res$BETA[,1], .05)
  ret[["beta_intercept_95q"]] <- quantile(res$BETA[,1], .95)
  for (i in 1:pr) {
    ret[[paste0("beta_row",i,"_mean")]] <- mean(res$BETA[,1+i])
    ret[[paste0("beta_row",i,"_5q")]] <- quantile(res$BETA[,1+i], .05)
    ret[[paste0("beta_row",i,"_95q")]] <- quantile(res$BETA[,1+i], .95)
  }
  for (j in 1:pc) {
    ret[[paste0("beta_col",j,"_mean")]] <- mean(res$BETA[,j+pr+1])
    ret[[paste0("beta_col",j,"_5q")]] <- quantile(res$BETA[,j+pr+1], .05)
    ret[[paste0("beta_col",j,"_95q")]] <- quantile(res$BETA[,j+pr+1], .95)
  }
  # 2. Projected estimates and CI's
  ret[["delta_intercept_mean"]] <- mean(res$DELTA[,1])
  ret[["delta_intercept_5q"]] <- quantile(res$DELTA[,1], .05)
  ret[["delta_intercept_95q"]] <- quantile(res$DELTA[,1], .95)
  for (i in 1:pr) {
    ret[[paste0("delta_row",i,"_mean")]] <- mean(res$DELTA[,1+i])
    ret[[paste0("delta_row",i,"_5q")]] <- quantile(res$DELTA[,1+i], .05)
    ret[[paste0("delta_row",i,"_95q")]] <- quantile(res$DELTA[,1+i], .95)
  }
  for (j in 1:pc) {
    ret[[paste0("delta_col",j,"_mean")]] <- mean(res$DELTA[,j+pr+1])
    ret[[paste0("delta_col",j,"_5q")]] <- quantile(res$DELTA[,j+pr+1], .05)
    ret[[paste0("delta_col",j,"_95q")]] <- quantile(res$DELTA[,j+pr+1], .95)
  }
  # 2. Difference in random effect estimates
  ret[["relmove_a_proj"]] <- sum((res$APM - res$APM.ORTH)^2)/sum(res$APM^2)
  ret[["relmove_b_proj"]] <- sum((res$BPM - res$BPM.ORTH)^2)/sum(res$BPM^2)
  # 3. Unprojected variance component estimates and CI's
  ret[["sigma_a_mean"]] <- mean(res$VC[,1])
  ret[["sigma_a_5q"]] <- quantile(res$VC[,1], .05)
  ret[["sigma_a_95q"]] <- quantile(res$VC[,1], .95)
  ret[["sigma_b_mean"]] <- mean(res$VC[,3])
  ret[["sigma_b_5q"]] <- quantile(res$VC[,3], .05)
  ret[["sigma_b_95q"]] <- quantile(res$VC[,3], .95)
  # 4. Projected variance component estimates and CI's
  ret[["sigma_a_proj_mean"]] <- mean(res$VC.ORTH[,1])
  ret[["sigma_a_proj_5q"]] <- quantile(res$VC.ORTH[,1], .05)
  ret[["sigma_a_proj_95q"]] <- quantile(res$VC.ORTH[,1], .95)
  ret[["sigma_b_proj_mean"]] <- mean(res$VC.ORTH[,3])
  ret[["sigma_b_proj_5q"]] <- quantile(res$VC.ORTH[,3], .05)
  ret[["sigma_b_proj_95q"]] <- quantile(res$VC.ORTH[,3], .95)
  
  # This is the output of the block:
  ret
}
output <- as.data.frame(output)

# Output the file
write.csv(output, paste0(outfilebase, run, ".csv"))

