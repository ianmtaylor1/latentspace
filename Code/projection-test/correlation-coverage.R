library(foreach)

n <- 20
beta <- 1
gamma <- 1

runs <- 10
cors <- c(-0.9, -0.5,  0, 0.5, 0.9)
reps <- 100

simulation.data <- foreach(i = 1:runs, .combine="rbind") %do% {
  cat("**** RUN = ", i, "\n", sep="")
  # Generate the hidden covariate (row covariate)
  X.unobserved <- rnorm(n)
  foreach(rho = cors, .combine="rbind") %do% {
    cat("** rho = ", rho, "\n", sep="")
    # Correlated predictor
    X <- c(amenproj::rmvnorm(1, rho*X.unobserved, (1-rho^2)*diag(n)))
    obs.rho <- cor(X, X.unobserved)
    foreach(j = 1:reps, .combine="rbind") %do% {
      cat("Rep = ", j, "\n", sep="")
      # Observations (Y)
      Z <- (X*beta + X.unobserved*gamma) %*% t(rep(1,n)) + matrix(rnorm(n^2), nrow=n, ncol=n)
      Y <- (Z >= 0) * 1.0
      # Fit it without projections
      fit <- amenproj::ame(Y=Y, Xrow=X, intercept=FALSE, family="bin", 
                           project=FALSE,
                           rvar=TRUE, cvar=FALSE, R=0, dcor=FALSE, 
                           plot=FALSE, gof=FALSE, print=FALSE)
      # Fit it with projections
      fitproj <- amenproj::ame(Y=Y, Xrow=X, intercept=FALSE, family="bin", 
                               project=TRUE,
                               rvar=TRUE, cvar=FALSE, R=0, dcor=FALSE, 
                               plot=FALSE, gof=FALSE, print=FALSE)
      # Output
      data.frame(run=i, rep=j, rho=rho, obs.rho=obs.rho, 
                 beta=beta, gamma=gamma, 
                 pctile=mean(fit$BETA <= beta), pctile.proj=mean(fitproj$BETA <= beta))
    }
  }
}