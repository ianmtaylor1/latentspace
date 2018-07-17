library(doParallel)

# Calculate the density of a MVN distribution
dmvnorm <- function(x,
                    mu = rep(0,length(x)),
                    Sigma = diag(length(x)),
                    Sigma.inv = chol2inv(chol(Sigma))
                   ) {
  sqrt(det(Sigma.inv) / (2 * pi)) * exp(-0.5 * t(x - mu) %*% Sigma.inv %*% (x - mu))
}

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
  
  # Arrays to store densities in
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


# Computes the bayes factor for comparing a model with no random effects vs
# a model with a row random effect.
# Only works for normal models, without an intercept, non-symmetric
# Does the same thing as the other function, but where it can, will use random
# Effects as part of the covariance of Y instead of having to sample over them.
# Maybe faster?
bayes_factor_covar <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, intercept=FALSE,
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
  
  # Convenient pieces to use later
  cY <- c(Y)                 # Columnized version of Y.
  cY_na_idx <- !is.na(c(Y))  # Boolean array indexing NA's in cy
  nn <- sum(cY_na_idx)       # Number of non-NA entries in Y
  
  # Build nn-by-nn matrices to use as a base for adding in covariance
  # structures into the covariance matrix for calculating normal density
  # as necessary
  # Note: These are meant to be the covariance matrix for the columnized version
  # of Y, c(Y). R does "column first" conversion.
  
  # rcor_base is a base for adding row correlation (i.e. from 'a' random effects)
  # It has a 1 in every entry for two edges (i,j) and (i,k) in the same row
  rcor_base <- matrix(0, nrow=n*n, ncol=n*n)
  # ccor_base is a base for adding column correlation (i.e. from 'b' random effects)
  # It has a 1 in every entry for two edges (i,j) and (k,j) that are in the same column
  ccor_base <- matrix(0,nrow=n*n, ncol=n*n)
  # tcor_base is a base for adding "transitive" correlations. i.e. due to sigma-squared-ab,
  # when one edge's receiver is another edge's sender.
  tcor_base <- matrix(0,nrow=n*n, ncol=n*n)
  # Fill the matrices
  for (i in 1:n) {
    # Construct a vector for row i
    onecol <- rep(0,n)
    onecol[i] <- 1
    inrow <- rep(onecol, n)
    # Construct a vector for column i
    incol <- c(rep(0, n*(i-1)), rep(1, n), rep(0, n*(n-i)))
    # Apply outer products and add to correlation base matrices
    rcor_base <- rcor_base + outer(inrow, inrow, "*")
    ccor_base <- ccor_base + outer(incol, incol, "*")
    tcor_base <- tcor_base + outer(inrow, incol, "*") + outer(incol, inrow, "*")
  }
  rcor_base <- rcor_base[cY_na_idx,cY_na_idx]
  ccor_base <- ccor_base[cY_na_idx,cY_na_idx]
  tcor_base <- tcor_base[cY_na_idx,cY_na_idx] # + 2 * diag(nn), if the amen vignette is to be believed
  
  # dcor_base is a base for adding dyad correlation. It has a 1 in every entry
  # for a dyad pair e.g. edges (i,j) and (j,i), and a zero everywhere else
  dcor_base <- matrix(0, nrow=n*n, ncol=n*n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dcor_base[n*(i-1)+j,n*(j-1)+i] <- 1
      dcor_base[n*(j-1)+9,n*(i-1)+j] <- 1
    }
  }
  dcor_base <- dcor_base[cY_na_idx,cY_na_idx]
  
  # Arrays to store densities in
  densities_no_rvar <- c()
  densities_rvar <- c()
  
  for (s in 1:samples) {
    # Draw sigma-squared from prior
    sigma2 <- 1/rgamma(n=1, shape=1/2, scale=1/2)
    # Draw sigma-squared-a from prior
    sigma2a <- 1/rgamma(n=1, shape=1/2, scale=1/2)
    # Draw beta from prior
    if (p > 0) {
      beta <- c(rmvnorm(n=1, mu=rep(0,p), Sigma=n^2 * XX))  # Note: if p==0, this will produce numeric(0)
    } else {
      beta <- numeric(0)
    }
    Xb <- Xbeta(X, beta)
    # Calculate contional density of Y for every non-NA element in matrix
    # using the dmvnorm function
    d_rvar <- dmvnorm(cY[cY_na_idx],
                      mu = c(Xb)[cY_na_idx],
                      Sigma = sigma2 * diag(nn) + sigma2a * rcor_base)
    d_no_rvar <- dmvnorm(cY[cY_na_idx],
                         mu = c(Xb)[cY_na_idx],
                         Sigma = sigma2 * diag(nn))
    
    densities_rvar <- c(densities_rvar, d_rvar)
    densities_no_rvar <- c(densities_no_rvar, d_no_rvar)
  }
  
  mean(densities_rvar)/mean(densities_no_rvar)
}

##### Main section of code

n <- 20
a_var <- 1
e_var <- 1
runs <- 20
reps <- 100
workers <- 2

rndseedbase <- 10000000

registerDoParallel(cl=workers, cores=workers)

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
    #data.frame(rep=rep, run=run,
    #           bf_no_rvar=bayes_factor_covar(y_no_rvar),
    #           bf_rvar=bayes_factor_covar(y_rvar))
  }
}

stopImplicitCluster()

write.csv(results,file="bayes_factor_results_covarversion.csv",row.names=FALSE)