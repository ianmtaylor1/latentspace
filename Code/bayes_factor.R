library(doParallel)

# Computes the bayes factor for comparing a model with no random effects vs
# a model with a row random effect.
# Only works for normal models, without an intercept, non-symmetric
bayes_factor <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, intercept=FALSE,
                         seed=1, samples=10000) {
  library(amen)
  # Baked in parameters for this version
  model <- "nrm"
  intercept <- FALSE
  prior <- list()
  
  set.seed(seed)
  
  # Build the design matrix
  diag(Y) <- NA
  n <- nrow(Y)
  X <- design_array(Xrow,Xcol,Xdyad,intercept,n)
  p <- dim(X)[3]
  
  XX <- t(apply(X,3,c)) %*% apply(X,3,c)

  densities_no_rvar <- c()
  densities_rvar <- c()
  
  for (s in 1:samples) {
    # Draw sigma-squared from prior
    sigma2 <- 1/rgamma(n=1, shape=1/2, scale=1/2)
    # Draw sigma-squared-a from prior
    sigma2a <- 1/rgamma(n=1, shape=1/2, scale=1/2)
    # Draw a from conditional prior, given sigma-squared-a
    a <- rnorm(n=n, mean=0, sd=sqrt(sigma2a))
    # Draw beta from prior
    if (p > 0) {
      beta <- c(rmvnorm(n=1, mu=rep(0,p), Sigma=n^2 * XX))  # Note: if p==0, this will produce numeric(0)
    } else {
      beta <- numeric(0)
    }
    Xb <- Xbeta(X, beta)
    # Calculate condition density of Y for every non-NA element in matrix
    # (is there a faster way to do this?)
    d_rvar <- 1
    d_no_rvar <- 1
    for (i in 1:n) {
      for (j in 1:n) {
        if (!is.na(Y[i,j])) {
          # Calculate density with a
          d_rvar <- d_rvar * dnorm(Y[i,j], mean=Xb[i,j]+a[i], sd=sqrt(sigma2))
          # Calculate density without a
          d_no_rvar <- d_no_rvar * dnorm(Y[i,j], mean=Xb[i,j], sd=sqrt(sigma2))
        }
      }
    }
    densities_rvar <- c(densities_rvar, d_rvar)
    densities_no_rvar <- c(densities_no_rvar, d_no_rvar)
  }
  
  mean(densities_rvar)/mean(densities_no_rvar)
}


##### Main section of code

n <- 20
a_var <- 1
e_var <- 1
runs <- 10
reps <- 100

rndseedbase <- 10000000

registerDoParallel(cl=4, cores=4)

results <- foreach(run=1:runs, .combine="rbind") %do% {
  cat("Run",run,"-",date(),"\n")
  set.seed(rndseedbase*run)
  # Generate a
  a <- rnorm(n, 0, a_var)
  # The output of this foreach gets returned since it's the last thing in this loop
  foreach(rep=1:reps, .combine="rbind") %dopar% {
    cat("Run",run,", Rep",rep,"-",date(),"\n")
    set.seed(rndseedbase*run + rep)
    # Generate errors
    e <- matrix(rnorm(n*n, 0, e_var), nrow=n, ncol=n)
    
    # Create two y's: with and without sender effects
    y_no_rvar <- e
    y_rvar <- e + t(sapply(a,function(x){rep(x,n)})) #Converts 'a' into a matrix with each element repeated for one row
    
    # Calculate the bayes factor and output it
    data.frame(rep=rep, run=run,
               bf_no_rvar=bayes_factor(y_no_rvar),
               bf_rvar=bayes_factor(y_rvar))
  }
}

stopImplicitCluster()

write.csv(results,file="bayes_factor_results.csv",row.names=FALSE)