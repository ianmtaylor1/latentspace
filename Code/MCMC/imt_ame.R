library(amen)

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

update_beta_fc <- function(Y, X, sigma.squared, a, b, sigma.ab) {
  # Hyperparameters
  n <- dim(Y)[1]
  p <- dim(X)[3]
  Xmat <- apply(X, 3, c)[nafilter(n),]
  XtX <- t(Xmat) %*% Xmat
  Sigma_beta <- solve(XtX) * (n ^ 2)
  # Base matrix used in both mean and covar of beta
  V <- solve(XtX + solve(Sigma_beta / sigma.squared))
  # Y without all other effects
  yprime <- Y - (a %*% t(rep(1, n))) - (rep(1, n) %*% t(b))
  ycol <- c(yprime)[nafilter(n)]
  # The mean of beta's full conditional distribution
  fcmean <- V %*% t(Xmat) %*% ycol
  # Generate new beta
  beta <- c(rmvnorm(n=1, mu=fcmean, Sigma=sigma.squared * V))
  return(beta)
}

update_sigma.squared_fc <- function(Y, X, beta, a, b, sigma.ab) {
  # Hyperparameters
  alpha_sigma <- 0.5
  beta_sigma <- 0.5
  # Find n
  n <- dim(Y)[1]
  # Y without all other effects
  yprime <- Y - (a %*% t(rep(1, n))) - (rep(1, n) %*% t(b)) - Xbeta(X, beta)
  ycol <- c(yprime)[nafilter(n)]
  # Parameters of fc gamma distribution of sigma.squared inverse
  shape <- n/2 + alpha_sigma
  rate <- (t(ycol) %*% ycol) / 2 + 1/beta_sigma
  # Generate new sigma.squared
  sigma.squared <- 1/rgamma(n=1, shape=shape, rate=rate)
  return(sigma.squared)
}

update_a_fc <- function(Y, X, beta, sigma.squared, b, sigma.ab) {
  # Calculate rowmeans of Y, without other effects
  yprime <- Y - (rep(1, n) %*% t(b)) - Xbeta(X, beta)
  n <- dim(Y)[1]
  yrowmeans <- rowSums(yprime, na.rm=TRUE) / (n - 1)
  # Calculate variance and meanvector of 
  var.a <- 1 / (((n - 1) / sigma.squared) + (sigma.ab[2,2] / det(sigma.ab)))
  mean.a <- var.a * (((n-1) * yrowmeans / sigma.squared) + (sigma.ab[1,2] * b / det(sigma.ab)))
  # Calculate new a's
  a <- rnorm(n=n, mean=mean.a, sd=sqrt(var.a))
  return(a)
}

update_b_fc <- function(Y, X, beta, sigma.squared, a, sigma.ab) {
  # Calculate colmeans of Y, without other effects
  yprime <- Y - (a %*% t(rep(1, n))) - Xbeta(X, beta)
  n <- dim(Y)[1]
  ycolmeans <- colSums(yprime, na.rm=TRUE) / (n - 1)
  # Calculate variance and meanvector of 
  var.b <- 1 / (((n - 1) / sigma.squared) + (sigma.ab[1,1] / det(sigma.ab)))
  mean.b <- var.b * (((n-1) * ycolmeans / sigma.squared) + (sigma.ab[1,2] * a / det(sigma.ab)))
  # Calculate new b's
  b <- rnorm(n=n, mean=mean.b, sd=sqrt(var.b))
  return(b)
}

update_sigma.ab_fc <- function(Y, X, beta, sigma.squared, a, b) {
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

#########################################################################
###### MAIN MCMC FUNCTION ###############################################
#########################################################################

imt_ame <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, intercept=TRUE,
                    burn=500, n.iter=10000, save.every=25) {
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
  
  # Lists to save the samples in temporarily
  samples.a <- list()
  samples.b <- list()
  samples.beta <- list()
  samples.sigma.squared <- list()
  samples.sigma.ab <- list()
  
  # Iterate through the gibbs sampler
  for (i in 1:(n.iter + burn)) {
    # Update sigma.squared
    sigma.squared <- update_sigma.squared_fc(Y=Y, X=X, beta=beta, a=a, b=b, sigma.ab=sigma.ab)
    # Update Sigma.ab
    sigma.ab <- update_sigma.ab_fc(Y=Y, X=X, beta=beta, sigma.squared=sigma.squared, a=a, b=b)
    # Update beta
    beta <- update_beta_fc(Y=Y, X=X, sigma.squared=sigma.squared, a=a, b=b, sigma.ab=sigma.ab)
    # Update a,b
    a <- update_a_fc(Y=Y, X=X, beta=beta, sigma.squared=sigma.squared, b=b, sigma.ab=sigma.ab)
    b <- update_b_fc(Y=Y, X=X, beta=beta, sigma.squared=sigma.squared, a=a, sigma.ab=sigma.ab)
    # Update latent y? (TBD: probit model)
    if ((i > burn) && ((i-burn) %% save.every == 0)) {
      cat("Saving...",i,"\n")
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





