
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
reps <- 50       # How many times to generate Y and fit?
pr <- 2          # How many fixed row covariates?
#pc <- 1          # How many fixed column covariates?
pd <- 2           # How many dyadic covariates?
#noise <- 1       # Noise level (magnitude of unobserved effects relative to cumulative effect of fixed)
#set.corr <- 0.5  # Canonical correlation of unobserved effects with fixed effects
n <- 27          # Number of nodes in networks
outfilebase <- NULL # Base of name of CSV to write to.
# Read the command line. May override variables set in defaults.
parseCommandArgs()

################################################################################

param.idx <- ((run - 1) %/% 4) + 1

# Parameters to use for each run
corr.all <- c(0.1, 0.1, 0.9, 0.9)
noise.all <- c(1, 2, 1, 2)
startreps.all <- c(1, 51, 101, 151)

# Choose based on run
noise <- noise.all[(param.idx - 1) %% 4 + 1]
set.corr <- corr.all[(param.idx-1) %% 4 + 1]
startrep <- startreps.all[(run - 1) %% 4 + 1]

# Repeatability
base.seed <- 111

set.seed(base.seed * 1000 + param.idx)

# Parameters for generating data
intercept <- -1
beta.r <- rep(1, pr)
beta.d <- rep(1, pd)

# Generating the fixed things for this run
Xr <- matrix(rnorm(pr*n), ncol=pr)
Xd <- array(rnorm(pd*n*n), dim=c(n,n,pd))
a.true <- gencancor(cbind(rep(1,n), Xr), rnorm(n), set.corr) * noise * sqrt(sum(beta.r^2))
b.true <- rep(0, n)

# MCMC controlling parameters
iter<-80000
burn<-5000
thin<-40

# What is the true delta for this set of parameters?
Xmat <- matrix(0, nrow=n^2, ncol=1+pd+pr)
Xmat[,1] <- 1
for (i in 1:pd) {
  Xmat[,1+i] <- c(Xd[,,i])
}
for (j in 1:pr) {
  Xmat[,1+pd+j] <- rep(Xr[,j], n)
}
a.full <- rep(a.true, n)
delta <- c(intercept, beta.d, beta.r) + c(solve(t(Xmat) %*% Xmat) %*% (t(Xmat) %*% a.full))
delta.int <- delta[1]
delta.d <- delta[2:(pd+1)]
delta.r <- delta[(pd+2):(1+pd+pr)]


output <- foreach(rep=startrep:(startrep+reps), .combine="rbind") %do% {
  cat("Rep",rep,"\n")
  
  set.seed((base.seed * 1000 + param.idx) * 10000 + rep)
  
  # Response variable
  dyadpart <- matrix(0, nrow=n, ncol=n)
  for (i in 1:pd) {
    dyadpart <- dyadpart + beta.d[i] * Xd[,,i]
  }
  Y <- (
    intercept
    + dyadpart
    + outer(c(Xr %*% beta.r), rep(0, n), "+")
    + outer(a.true, b.true, "+")
    + matrix(rnorm(n*n), nrow=n) 
    >= 0)
  
  # Run the MCMC
  res <- ame(Y, Xrow=Xr, Xdyad=Xd, family="bin", halfcauchy=TRUE, project=TRUE, rvar=TRUE, cvar=FALSE, dcor=FALSE,
             nscan=iter, burn=burn, odens=thin,
             print=FALSE, plot=FALSE, gof=FALSE)
  res.no.re <- ame(Y, Xrow=Xr, Xdyad=Xd, family="bin", rvar=FALSE, cvar=FALSE, dcor=FALSE,
                   nscan=iter, burn=burn, odens=thin,
                   print=FALSE, plot=FALSE, gof=FALSE)
  
  # Analyze results
  ret <- data.frame(run=param.idx, rep=rep, pr=pr, pd=pd, noise=noise, corr=set.corr, n=n, stringsAsFactors=FALSE)
  # 0. True parameter values
  ret[["beta_intercept_true"]] <- intercept
  ret[["delta_intercept_true"]] <- delta.int
  for (i in 1:pr) {
    ret[[paste0("beta_row",i,"_true")]] <- beta.r[i]
    ret[[paste0("delta_row",i,"_true")]] <- delta.r[i]
  }
  for (i in 1:pd) {
    ret[[paste0("beta_dyad",i,"_true")]] <- beta.d[i]
    ret[[paste0("delta_dyad",i,"_true")]] <- delta.d[i]
  }
  # 1. Unprojected estimates and CI's
  ret[["beta_intercept_mean"]] <- mean(res$BETA[,1])
  ret[["beta_intercept_5q"]] <- quantile(res$BETA[,1], .05)
  ret[["beta_intercept_95q"]] <- quantile(res$BETA[,1], .95)
  for (i in 1:pr) {
    ret[[paste0("beta_row",i,"_mean")]] <- mean(res$BETA[,1+i])
    ret[[paste0("beta_row",i,"_5q")]] <- quantile(res$BETA[,1+i], .05)
    ret[[paste0("beta_row",i,"_95q")]] <- quantile(res$BETA[,1+i], .95)
  }
  for (j in 1:pd) {
    ret[[paste0("beta_dyad",j,"_mean")]] <- mean(res$BETA[,j+pr+1])
    ret[[paste0("beta_dyad",j,"_5q")]] <- quantile(res$BETA[,j+pr+1], .05)
    ret[[paste0("beta_dyad",j,"_95q")]] <- quantile(res$BETA[,j+pr+1], .95)
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
  for (j in 1:pd) {
    ret[[paste0("delta_dyad",j,"_mean")]] <- mean(res$DELTA[,j+pr+1])
    ret[[paste0("delta_dyad",j,"_5q")]] <- quantile(res$DELTA[,j+pr+1], .05)
    ret[[paste0("delta_dyad",j,"_95q")]] <- quantile(res$DELTA[,j+pr+1], .95)
  }
  # 3. Beta estimates with no random effect
  ret[["no_re_beta_intercept_mean"]] <- mean(res.no.re$BETA[,1])
  ret[["no_re_beta_intercept_5q"]] <- quantile(res.no.re$BETA[,1], .05)
  ret[["no_re_beta_intercept_95q"]] <- quantile(res.no.re$BETA[,1], .95)
  for (i in 1:pr) {
    ret[[paste0("no_re_beta_row",i,"_mean")]] <- mean(res.no.re$BETA[,1+i])
    ret[[paste0("no_re_beta_row",i,"_5q")]] <- quantile(res.no.re$BETA[,1+i], .05)
    ret[[paste0("no_re_beta_row",i,"_95q")]] <- quantile(res.no.re$BETA[,1+i], .95)
  }
  for (j in 1:pd) {
    ret[[paste0("no_re_beta_dyad",j,"_mean")]] <- mean(res.no.re$BETA[,j+pr+1])
    ret[[paste0("no_re_beta_dyad",j,"_5q")]] <- quantile(res.no.re$BETA[,j+pr+1], .05)
    ret[[paste0("no_re_beta_dyad",j,"_95q")]] <- quantile(res.no.re$BETA[,j+pr+1], .95)
  }
  # 4. Unprojected variance component estimates and CI's
  ret[["sigma_a_mean"]] <- mean(res$VC[,1])
  ret[["sigma_a_5q"]] <- quantile(res$VC[,1], .05)
  ret[["sigma_a_95q"]] <- quantile(res$VC[,1], .95)
  # 5. Projected variance component estimates and CI's
  ret[["sigma_a_proj_mean"]] <- mean(res$VC.ORTH[,1])
  ret[["sigma_a_proj_5q"]] <- quantile(res$VC.ORTH[,1], .05)
  ret[["sigma_a_proj_95q"]] <- quantile(res$VC.ORTH[,1], .95)
  
  # This is the output of the block:
  ret
}
output <- as.data.frame(output)

# Output the file
write.csv(output, paste0(outfilebase, run, ".csv"))

