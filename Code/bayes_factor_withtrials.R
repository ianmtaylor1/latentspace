library(doParallel)

# Calculate the density of a MVN distribution
dmvnorm <- function(x,
                    mu = rep(0,length(x)),
                    Sigma = diag(length(x))
                   ) {
  n <- length(x)
  Sigma.chol <- chol(Sigma)
  sqrt.det.Sigma <- prod(diag(Sigma.chol))
  inv.sqrt.2pi <- 0.39894228040143267794
  ((inv.sqrt.2pi ^ n) / sqrt.det.Sigma) * exp(-0.5 * t(x - mu) %*% chol2inv(Sigma.chol) %*% (x - mu))
}

# Computes the bayes factor for comparing a model with no random effects vs
# a model with a row random effect.
# Only works for normal models, without an intercept, non-symmetric
bayes_factor <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, intercept=FALSE,
                         prior_rvar=0.5, prior_cvar=0.5,
                         seed=NULL, samples=10000) {
  library(amen)
  # Baked in parameters for this version
  model <- "nrm"
  intercept <- FALSE
  prior <- list()
  
  if (!is.null(seed)) set.seed(seed)
  
  # Build the design matrix
  diag(Y) <- NA
  n <- nrow(Y)
  X <- design_array(Xrow,Xcol,Xdyad,intercept,n)
  p <- dim(X)[3]
  
  XX <- t(apply(X,3,c)) %*% apply(X,3,c)
  
  # Arrays to store densities in
  densities_none <- c()
  densities_a <- c()
  densities_b <- c()
  densities_ab <- c()
  
  for (s in 1:samples) {
    # Draw sigma-squared from prior
    sigma2 <- 1/rgamma(n=1, shape=1/2, scale=1/2)
    # Draw sigma-squared-a from prior (for only row effects)
    sigma2a <- 1/rgamma(n=1, shape=1/2, scale=1/2)
    # Draw sigma-squared-b from prior (for only column effects)
    sigma2b <- 1/rgamma(n=1, shape=1/2, scale=1/2)
    # Draw a from conditional prior, given sigma-squared-a
    a <- rnorm(n=n, mean=0, sd=sqrt(sigma2a))
    # Draw b from conditional prior, given sigma-squared-b
    b <- rnorm(n=n, mean=0, sd=sqrt(sigma2b))
    # Draw sigma-squared-ab from prior (for both row and column effects)
    sigma2ab <- solve(rwish(S0=diag(2), nu=3))
    # Draw ab from conditional prior, given sigma-squared-ab
    ab <- rmvnorm(n=n, mu=rep(0,2), Sigma=sigma2ab)
    # Draw beta from prior
    if (p > 0) {
      beta <- c(rmvnorm(n=1, mu=rep(0,p), Sigma=n^2 * XX))  # Note: if p==0, this will produce numeric(0)
    } else {
      beta <- numeric(0)
    }
    Xb <- Xbeta(X, beta)
    # Calculate condition density of Y for every non-NA element in matrix
    # (is there a faster way to do this?)
    d_none <- 1
    d_a <- 1
    d_b <- 1
    d_ab <- 1
    for (i in 1:n) {
      for (j in 1:n) {
        if (!is.na(Y[i,j])) {
          # Calculate density for no row or column effects
          d_none <- d_none * dnorm(Y[i,j], mean=Xb[i,j], sd=sqrt(sigma2))
          # Calculate density with a row effect only
          d_a <- d_a * dnorm(Y[i,j], mean=Xb[i,j]+a[i], sd=sqrt(sigma2))
          # Calculate density with a column effect only
          d_b <- d_b * dnorm(Y[i,j], mean=Xb[i,j]+b[j], sd=sqrt(sigma2))
          # Calculate density with both a row and column effect
          d_ab <- d_ab * dnorm(Y[i,j], mean=Xb[i,j]+ab[i,1]+ab[j,2], sd=sqrt(sigma2))
        }
      }
    }
    densities_none <- c(densities_none, d_none)
    densities_a <- c(densities_a, d_a)
    densities_b <- c(densities_b, d_b)
    densities_ab <- c(densities_ab, d_ab)
  }
  
  avg_none <- mean(densities_none)
  avg_a <- mean(densities_a)
  avg_b <- mean(densities_b)
  avg_ab <- mean(densities_ab)
  
  totalweight <- prior_rvar * prior_cvar * avg_ab +
                 (1 - prior_rvar) * prior_cvar * avg_b +
                 prior_rvar * (1 - prior_cvar) * avg_a +
                 (1 - prior_rvar) * (1 - prior_cvar) * avg_none
  
  posterior <- c(
    (1 - prior_rvar) * (1 - prior_cvar) * avg_none / totalweight,
    prior_rvar * (1 - prior_cvar) * avg_a / totalweight,
    (1 - prior_rvar) * prior_cvar * avg_b / totalweight,
    prior_rvar * prior_cvar * avg_ab / totalweight
  )
  names(posterior) <- c("none","rvar","cvar","rcvar")
  
  posterior
}

##### Main section of code

n <- 20
a_var <- 1
e_var <- 1
runs <- c(7,10)
reps <- 1:100
trials <- 1:20
workers <- 2

rndseedbase <- 10000000

registerDoParallel(cl=workers, cores=workers)

results <- foreach(run=runs, .combine="rbind") %do% {
  cat("Run",run,"-",date(),"\n")
  set.seed(rndseedbase*run)
  # Generate a
  a <- rnorm(n, 0, a_var)
  # The output of this foreach gets returned since it's the last thing in this loop
  foreach(rep=reps, .combine="rbind") %do% {
    cat("Run",run,", Rep",rep,"-",date(),"\n")
    set.seed(rndseedbase*run + rep)
    # Generate errors
    e <- matrix(rnorm(n*n, 0, e_var), nrow=n, ncol=n)
    
    # Create two y's: with and without sender effects
    y_no_rvar <- e
    y_rvar <- e + t(sapply(a,function(x){rep(x,n)})) #Converts 'a' into a matrix with each element repeated for one row
    
    # Repeat each "rep" for a few "trials", with a new seed at the start of the bayes_factor
    # routine. This is to test consistency/stability
    foreach(trial=trials, .combine="rbind") %dopar% {
      # Calculate the bayes factor and output it
      bf_no_rvar <- bayes_factor(y_no_rvar, seed=trial)
      bf_rvar <- bayes_factor(y_rvar, seed=trial)
      data.frame(rep=rep, run=run, trial=trial,
                 bf_no_rvar=bf_no_rvar["rvar"]/bf_no_rvar["none"],
                 bf_rvar=bf_rvar["rvar"]/bf_rvar["none"])
    }
  }
}

stopImplicitCluster()

write.csv(results,file="bayes_factor_results_consistency.csv",row.names=FALSE)