
library(amenhs)

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

# Repeatability
base.seed <- 12345
set.seed(base.seed)

# Parameters for generating data
n <- 27
pr <- 5
pc <- 1
intercept <- -1
beta.r <- rep(1, pr)
beta.c <- rep(1, pc)

noise <- 1 # Approximate Noise-to-signal ratio
set.corr <- 0.9 # Desired canonical correlation of unobserved effects to observed

# Generating the fixed things for this run
Xr <- matrix(rnorm(pr*n), ncol=pr)
Xc <- matrix(rnorm(pc*n), ncol=pc)
a.true <- gencancor(cbind(rep(1,n), Xr), rnorm(n), set.corr) * noise * sqrt(sum(beta.r^2))
b.true <- gencancor(cbind(rep(1,n), Xc), rnorm(n), set.corr) * noise * sqrt(sum(beta.c^2))

# MCMC controlling parameters
iter<-100000
burn<-5000
thin<-50

# How many reps to generate Y for
reps <- 20

for (rep in 1:reps) {
  set.seed((reps+1)*base.seed + rep)
  # Response variable
  Y <- (intercept + outer(c(Xr %*% beta.r), c(Xc %*% beta.c), "+") + outer(a.true, b.true, "+") + matrix(rnorm(n*n), nrow=n) >= 0)
  
  # Run the MCMC
  res <- ame(Y, Xrow=Xr, Xcol=Xc, family="bin", halfcauchy=TRUE, project=TRUE, dcor=FALSE,
             nscan=iter, burn=burn, odens=thin,
             print=TRUE, plot=FALSE, gof=FALSE)
  
  # Analyze results
  colMeans(res$BETA)
  apply(res$BETA, MARGIN=2, FUN=quantile, probs=c(0.05, 0.95))
  colMeans(res$DELTA)
  apply(res$DELTA, MARGIN=2, FUN=quantile, probs=c(0.05, 0.95))
  # 1. Difference in beta estimates
  # 2. Difference in random effect estimates
  # 3. Difference in random effect variance component estimates
}

